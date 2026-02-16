#!/usr/bin/env python3
"""
Stan Variable Usage Analyzer

Traces variable usage through Stan model files to identify:
1. Declared but unused variables
2. Variables used only in transformed data
3. Variables used only in generated quantities
"""

import re
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Tuple
import sys


class StanVariableAnalyzer:
    def __init__(self, stan_dir: Path):
        self.stan_dir = Path(stan_dir)
        self.variables = defaultdict(lambda: {
            'declared_in': [],
            'used_in': [],
            'usage_contexts': []
        })

        # Track include chains
        self.include_cache = {}

    def resolve_includes(self, file_path: Path, base_dir: Path = None) -> List[Path]:
        """Recursively resolve all #include directives in a Stan file."""
        if base_dir is None:
            base_dir = file_path.parent

        file_path = file_path.resolve()
        if file_path in self.include_cache:
            return self.include_cache[file_path]

        resolved = [file_path]

        try:
            content = file_path.read_text()
        except Exception as e:
            print(f"Warning: Could not read {file_path}: {e}")
            return resolved

        # Find all #include directives
        include_pattern = r'#include\s+"([^"]+)"'
        for match in re.finditer(include_pattern, content):
            include_path = match.group(1)

            # Resolve relative to base_dir
            full_path = (base_dir / include_path).resolve()

            if full_path.exists() and full_path not in resolved:
                # Recursively resolve includes in the included file
                resolved.extend(self.resolve_includes(full_path, base_dir))

        self.include_cache[file_path] = resolved
        return resolved

    def extract_data_variables(self, content: str, context: str) -> Set[str]:
        """Extract variable declarations from data block."""
        variables = set()

        # Match various Stan data declaration patterns
        patterns = [
            r'(?:int|real|vector|matrix|array\[.*?\]\s+(?:int|real))\s*(?:<[^>]+>)?\s+(\w+)',
            r'array\[.*?\]\s+(\w+)\s*<',  # array[N] foo<...>
        ]

        for pattern in patterns:
            for match in re.finditer(pattern, content):
                var_name = match.group(1)
                if var_name and not var_name.startswith('_'):  # Skip internal vars
                    variables.add(var_name)
                    self.variables[var_name]['declared_in'].append(context)

        return variables

    def extract_transformed_data_variables(self, content: str, context: str) -> Set[str]:
        """Extract variable declarations/assignments from transformed data block."""
        variables = set()

        # Match various Stan transformed data declaration/assignment patterns
        patterns = [
            # Type declarations: int x = ...; or int<...> x = ...;
            r'^\s*(?:int|real|vector|matrix|row_vector|array\[.*?\]\s*(?:int|real|vector|matrix))\s*(?:<[^>]+>)?\s+(\w+)\s*[=;]',
            # Array declarations: array[...] type name
            r'^\s*array\[[^\]]+\]\s+(?:int|real|vector|matrix|row_vector)\s*(?:<[^>]+>)?\s+(\w+)\s*[=;]',
        ]

        lines = content.split('\n')
        for line_num, line in enumerate(lines, 1):
            for pattern in patterns:
                for match in re.finditer(pattern, line):
                    var_name = match.group(1)
                    if var_name and not var_name.startswith('_'):
                        variables.add(var_name)
                        self.variables[var_name]['declared_in'].append(f"{context}:{line_num}")

        return variables

    def find_variable_usage(self, content: str, var_name: str, context: str) -> List[str]:
        """Find all usages of a variable in content."""
        usages = []

        # Look for variable usage (word boundary to avoid partial matches)
        pattern = r'\b' + re.escape(var_name) + r'\b'

        lines = content.split('\n')
        for line_num, line in enumerate(lines, 1):
            if re.search(pattern, line):
                # Skip if it's just the declaration line
                if not re.match(r'^\s*(?:int|real|vector|matrix|array)', line):
                    usages.append(f"{context}:{line_num}: {line.strip()[:80]}")

        return usages

    def categorize_context(self, file_path: Path) -> str:
        """Categorize what section of the model a file belongs to."""
        path_str = str(file_path)

        # Check transformed_data FIRST before checking _data
        if 'transformed_data.stan' in path_str:
            return 'transformed_data'
        elif 'transformed_parameters.stan' in path_str:
            return 'transformed_parameters'
        elif '_data.stan' in path_str or 'base_data.stan' in path_str:
            return 'data'
        elif 'hyperparams.stan' in path_str or 'flags.stan' in path_str:
            return 'data'
        elif 'parameters.stan' in path_str:
            return 'parameters'
        elif 'priors.stan' in path_str:
            return 'model'
        elif 'generated_quantities' in path_str:
            return 'generated_quantities'
        elif file_path.suffix == '.stan':
            # Try to determine from content structure
            try:
                content = file_path.read_text()
                if 'model {' in content:
                    return 'model'
                elif 'generated quantities {' in content:
                    return 'generated_quantities'
            except:
                pass

        return 'unknown'

    def analyze_model(self, model_file: Path, analyze_type='data'):
        """Analyze a complete Stan model and all its includes.

        Args:
            model_file: Path to the main Stan model file
            analyze_type: 'data' or 'transformed_data'
        """
        print(f"\n{'='*80}")
        print(f"Analyzing: {model_file.name}")
        if analyze_type == 'transformed_data':
            print("Mode: Finding unused TRANSFORMED DATA variables")
        else:
            print("Mode: Finding unused DATA variables")
        print(f"{'='*80}\n")

        # Reset for this model
        self.variables.clear()
        self.include_cache.clear()

        # Resolve all includes
        all_files = self.resolve_includes(model_file, self.stan_dir)
        print(f"Found {len(all_files)} files (including transitive includes)\n")

        # Step 1: Extract all declared variables
        target_vars = set()

        if analyze_type == 'transformed_data':
            # Extract variables created in transformed data blocks
            for file_path in all_files:
                context_type = self.categorize_context(file_path)
                if context_type == 'transformed_data':
                    try:
                        content = file_path.read_text()
                        try:
                            rel_path = str(file_path.relative_to(self.stan_dir.resolve()))
                        except ValueError:
                            rel_path = file_path.name

                        vars_in_file = self.extract_transformed_data_variables(content, rel_path)
                        target_vars.update(vars_in_file)
                    except Exception as e:
                        print(f"Warning: Could not process {file_path}: {e}")
            print(f"Found {len(target_vars)} transformed data variables\n")
        else:
            # Extract variables declared in data blocks (original behavior)
            for file_path in all_files:
                context_type = self.categorize_context(file_path)
                if context_type == 'data':
                    try:
                        content = file_path.read_text()
                        try:
                            rel_path = str(file_path.relative_to(self.stan_dir.resolve()))
                        except ValueError:
                            rel_path = file_path.name

                        vars_in_file = self.extract_data_variables(content, rel_path)
                        target_vars.update(vars_in_file)
                    except Exception as e:
                        print(f"Warning: Could not process {file_path}: {e}")
            print(f"Found {len(target_vars)} data variables\n")

        # Step 2: Find usage of each variable
        for var_name in target_vars:
            for file_path in all_files:
                context_type = self.categorize_context(file_path)

                # Skip the source files themselves
                if analyze_type == 'transformed_data':
                    # When analyzing transformed_data vars, skip transformed_data files
                    if context_type == 'transformed_data':
                        continue
                else:
                    # When analyzing data vars, skip data files
                    if context_type == 'data':
                        continue

                try:
                    content = file_path.read_text()
                    # Try to get relative path, fallback to name if fails
                    try:
                        rel_path = str(file_path.relative_to(self.stan_dir.resolve()))
                    except ValueError:
                        rel_path = file_path.name

                    usages = self.find_variable_usage(content, var_name, rel_path)

                    if usages:
                        self.variables[var_name]['used_in'].extend(usages)
                        self.variables[var_name]['usage_contexts'].append(context_type)

                except Exception as e:
                    print(f"Warning: Could not search {file_path}: {e}")

        # Step 3: Categorize and report
        unused = []
        transformed_data_only = []
        gen_quantities_only = []
        used_in_model = []

        for var_name in sorted(target_vars):
            var_info = self.variables[var_name]
            contexts = set(var_info['usage_contexts'])

            if not contexts:
                unused.append(var_name)
            elif contexts == {'transformed_data'} and analyze_type != 'transformed_data':
                transformed_data_only.append(var_name)
            elif contexts == {'generated_quantities'}:
                gen_quantities_only.append(var_name)
            elif 'model' in contexts or 'transformed_parameters' in contexts or 'parameters' in contexts:
                used_in_model.append(var_name)
            else:
                # Mixed usage
                used_in_model.append(var_name)

        # Print results
        print(f"\n{'='*80}")
        print(f"RESULTS FOR: {model_file.name}")
        print(f"{'='*80}\n")

        if unused:
            print(f"🔴 UNUSED VARIABLES ({len(unused)}):")
            print("   (Never referenced anywhere)")
            for var in unused:
                decl_files = ', '.join(self.variables[var]['declared_in'])
                print(f"   - {var:30s} [declared in: {decl_files}]")
            print()

        if transformed_data_only:
            print(f"🟡 USED ONLY IN TRANSFORMED DATA ({len(transformed_data_only)}):")
            print("   (May be dead code if output is unused)")
            for var in transformed_data_only:
                print(f"   - {var:30s}")
            print()

        if gen_quantities_only:
            print(f"🟢 USED ONLY IN GENERATED QUANTITIES ({len(gen_quantities_only)}):")
            print("   (Optional analysis infrastructure)")
            for var in gen_quantities_only:
                print(f"   - {var:30s}")
            print()

        print(f"✅ USED IN MODEL/PARAMETERS ({len(used_in_model)}):")
        print(f"   (Active variables affecting model estimation)\n")

        # Summary stats
        total = len(target_vars)
        if analyze_type == 'transformed_data':
            print(f"Summary: {len(unused)}/{total} never used, "
                  f"{len(gen_quantities_only)}/{total} GQ only, "
                  f"{len(used_in_model)}/{total} used in model/params")
        else:
            print(f"Summary: {len(unused)}/{total} unused, "
                  f"{len(transformed_data_only)}/{total} transformed_data only, "
                  f"{len(gen_quantities_only)}/{total} generated_quantities only, "
                  f"{len(used_in_model)}/{total} used in model")

        return {
            'unused': unused,
            'transformed_data_only': transformed_data_only,
            'gen_quantities_only': gen_quantities_only,
            'used_in_model': used_in_model
        }


def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_stan_variables.py <stan_directory> [--transformed-data]")
        print("\nExample: python analyze_stan_variables.py /mnt/code/stan")
        print("         python analyze_stan_variables.py /mnt/code/stan --transformed-data")
        sys.exit(1)

    stan_dir = Path(sys.argv[1])
    analyze_type = 'transformed_data' if '--transformed-data' in sys.argv else 'data'

    if not stan_dir.exists():
        print(f"Error: Directory {stan_dir} does not exist")
        sys.exit(1)

    analyzer = StanVariableAnalyzer(stan_dir)

    # Find main model files
    models = [
        stan_dir / 'sf-ssm-log-space.stan',
        stan_dir / 'sf-ssls-lfo.stan',
        stan_dir / 'sf-ssls-lfo-endpoints.stan'
    ]

    all_results = {}
    for model in models:
        if model.exists():
            results = analyzer.analyze_model(model, analyze_type=analyze_type)
            all_results[model.name] = results
        else:
            print(f"Warning: Model {model} not found")

    # Print cross-model comparison
    print(f"\n\n{'='*80}")
    print("CROSS-MODEL COMPARISON")
    print(f"{'='*80}\n")

    # Find variables unused in ALL models
    if all_results:
        unused_in_all = set(all_results[list(all_results.keys())[0]]['unused'])
        for model_name, results in all_results.items():
            unused_in_all &= set(results['unused'])

        if unused_in_all:
            print(f"Variables unused in ALL models ({len(unused_in_all)}):")
            for var in sorted(unused_in_all):
                print(f"   - {var}")
        else:
            print("No variables are unused in all models.")


if __name__ == '__main__':
    main()
