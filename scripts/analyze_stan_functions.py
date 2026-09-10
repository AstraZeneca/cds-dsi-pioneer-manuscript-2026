#!/usr/bin/env python3
"""
Analyze Stan function definitions and usage to identify dead functions.
"""

import re
from pathlib import Path
from typing import Dict, List, Set, Tuple
from collections import defaultdict

class StanFunctionAnalyzer:
    def __init__(self, stan_dir: Path):
        self.stan_dir = stan_dir
        self.functions_defined: Dict[str, List[Tuple[Path, int]]] = defaultdict(list)  # func_name -> [(file, line)]
        self.functions_called: Dict[str, List[Tuple[Path, int, str]]] = defaultdict(list)  # func_name -> [(file, line, context)]
        self.all_content = ""

    def resolve_includes(self, file_path: Path, base_dir: Path = None) -> List[Path]:
        """Recursively resolve all #include directives."""
        if base_dir is None:
            base_dir = file_path.parent

        resolved = []
        with open(file_path, 'r') as f:
            content = f.read()

        # Find all #include statements
        for match in re.finditer(r'^\s*#include\s+"([^"]+)"', content, re.MULTILINE):
            include_path = match.group(1)
            full_path = (base_dir / include_path).resolve()

            if full_path.exists():
                resolved.append(full_path)
                # Recursively resolve includes in the included file
                resolved.extend(self.resolve_includes(full_path, base_dir))

        return resolved

    def extract_function_definitions(self, file_path: Path) -> List[Tuple[str, int]]:
        """Extract function definitions from a Stan file."""
        definitions = []

        with open(file_path, 'r') as f:
            content = f.read()
            lines = content.split('\n')

        # Pattern for function definitions (simplified, may need refinement)
        # Matches: return_type function_name(params) { or return_type function_name(params);
        pattern = r'^\s*(?:real|int|vector|row_vector|matrix|array|void)\s+(\w+)\s*\('

        for i, line in enumerate(lines, 1):
            # Skip comments
            if line.strip().startswith('//'):
                continue

            match = re.match(pattern, line)
            if match:
                func_name = match.group(1)
                # Skip Stan built-in functions that might match the pattern
                if func_name not in ['print', 'reject', 'fatal_error', 'target', 'get_lp']:
                    definitions.append((func_name, i))

        return definitions

    def find_function_calls(self, file_path: Path, func_name: str) -> List[Tuple[int, str]]:
        """Find all calls to a specific function in a file."""
        calls = []

        with open(file_path, 'r') as f:
            content = f.read()
            lines = content.split('\n')

        # Pattern for function calls: func_name( or func_name, or func_name)
        # The latter two catch function references (e.g., passed to map_rect)
        pattern = r'\b' + re.escape(func_name) + r'(?:\s*\(|[,\)])'

        for i, line in enumerate(lines, 1):
            # Skip comments
            if '//' in line:
                line = line[:line.index('//')]

            if re.search(pattern, line):
                calls.append((i, line.strip()))

        return calls

    def analyze_all_files(self, model_files: List[Path]):
        """Analyze all Stan files to build function definition and usage maps."""

        # Get all Stan files (not just models), excluding legacy directory
        all_stan_files = set()
        for model_file in model_files:
            all_stan_files.add(model_file)
            includes = self.resolve_includes(model_file)
            # Filter out legacy files
            includes = [f for f in includes if 'legacy' not in f.parts]
            all_stan_files.update(includes)

        # Extract function definitions from function files
        function_files = [
            self.stan_dir / "util.stan",
            self.stan_dir / "pos.stan",
            self.stan_dir / "gp.stan",
            self.stan_dir / "pfs_functions.stan",
            self.stan_dir / "lfo.stan",
            self.stan_dir / "recist.stanfunctions",
            self.stan_dir / "_sf_functions.stan"
        ]

        for func_file in function_files:
            if func_file.exists():
                definitions = self.extract_function_definitions(func_file)
                for func_name, line_num in definitions:
                    self.functions_defined[func_name].append((func_file, line_num))

        # Find function calls across all files (excluding legacy)
        for stan_file in all_stan_files:
            if 'legacy' in stan_file.parts:
                continue

            for func_name in self.functions_defined.keys():
                calls = self.find_function_calls(stan_file, func_name)
                for line_num, context in calls:
                    # Don't count the definition itself as a call
                    is_definition = any(
                        def_file == stan_file and def_line == line_num
                        for def_file, def_line in self.functions_defined[func_name]
                    )
                    if not is_definition:
                        self.functions_called[func_name].append((stan_file, line_num, context))

    def find_dead_functions(self) -> Dict[Path, List[Tuple[str, int]]]:
        """Find functions that are defined but never called."""
        dead_functions = defaultdict(list)

        for func_name, definitions in self.functions_defined.items():
            if func_name not in self.functions_called or len(self.functions_called[func_name]) == 0:
                for file_path, line_num in definitions:
                    dead_functions[file_path].append((func_name, line_num))

        return dead_functions

    def print_report(self):
        """Print a report of dead functions."""
        dead_funcs = self.find_dead_functions()

        total_defined = sum(len(defs) for defs in self.functions_defined.values())
        total_dead = sum(len(funcs) for funcs in dead_funcs.values())

        print(f"\n{'='*80}")
        print(f"STAN DEAD FUNCTION ANALYSIS")
        print(f"{'='*80}\n")

        print(f"Total functions defined: {total_defined}")
        print(f"Total functions used: {total_defined - total_dead}")
        print(f"Total DEAD functions: {total_dead} ({100*total_dead/total_defined:.1f}%)\n")

        if dead_funcs:
            print(f"{'='*80}")
            print("DEAD FUNCTIONS BY FILE:")
            print(f"{'='*80}\n")

            for file_path in sorted(dead_funcs.keys()):
                funcs = dead_funcs[file_path]
                rel_path = file_path.relative_to(self.stan_dir.parent)

                all_funcs_in_file = [
                    func_name for func_name, defs in self.functions_defined.items()
                    if any(def_file == file_path for def_file, _ in defs)
                ]

                pct = 100 * len(funcs) / len(all_funcs_in_file) if all_funcs_in_file else 0

                print(f"{rel_path}")
                print(f"  Dead: {len(funcs)}/{len(all_funcs_in_file)} functions ({pct:.1f}%)\n")

                for func_name, line_num in sorted(funcs, key=lambda x: x[1]):
                    print(f"    Line {line_num:4d}: {func_name}()")

                print()
        else:
            print("No dead functions found! ✅\n")


def main():
    stan_dir = Path("/mnt/code/stan")

    # Main model files
    model_files = [
        stan_dir / "sf-ssm-log-space.stan",
        stan_dir / "sf-ssls-lfo.stan",
        stan_dir / "sf-ssls-lfo-endpoints.stan"
    ]

    analyzer = StanFunctionAnalyzer(stan_dir)
    analyzer.analyze_all_files(model_files)
    analyzer.print_report()

    # Print summary by file for moving to legacy
    dead_funcs = analyzer.find_dead_functions()
    if dead_funcs:
        print(f"\n{'='*80}")
        print("RECOMMENDED ACTION:")
        print(f"{'='*80}\n")
        print("Move the following function files to stan/legacy/:\n")

        for file_path in sorted(dead_funcs.keys()):
            funcs = dead_funcs[file_path]
            all_funcs_in_file = [
                func_name for func_name, defs in analyzer.functions_defined.items()
                if any(def_file == file_path for def_file, _ in defs)
            ]

            pct = 100 * len(funcs) / len(all_funcs_in_file) if all_funcs_in_file else 0

            if pct > 50:  # More than 50% dead
                rel_path = file_path.relative_to(stan_dir.parent)
                print(f"  - {rel_path} (>{pct:.0f}% dead functions)")


if __name__ == "__main__":
    main()
