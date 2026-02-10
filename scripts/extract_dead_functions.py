#!/usr/bin/env python3
"""
Extract dead functions from Stan files and prepare them for removal.
"""

import re
from pathlib import Path
from typing import Dict, List, Set, Tuple
from collections import defaultdict

class StanFunctionExtractor:
    def __init__(self, stan_dir: Path):
        self.stan_dir = stan_dir
        self.dead_functions = {
            'stan/_sf_functions.stan': [
                ('sf_log_space_trajectory_lpdf', 151),
                ('calc_patient_states_rect', 611),
                ('assert_matching_states', 925)
            ],
            'stan/gp.stan': [
                ('ncp_gp_matern32', 143),
                ('ncp_gp_matern32', 150),
                ('ncp_gp_matern32', 157),
                ('ncp_gp_matern52', 164),
                ('ncp_gp_matern52', 171),
                ('multi_normal_lcdf', 311)
            ],
            'stan/pfs_functions.stan': [
                ('calc_pfs_n', 564),
                ('calc_admin_brier_score', 925)
            ],
            'stan/pos.stan': [
                ('get_sub_vector', 372),
                ('get_sub_row_vector', 387),
                ('get_last_int', 497)
            ],
            'stan/util.stan': [
                ('num_leq', 291),
                ('assert_greater', 687)
            ]
        }

    def extract_function_with_docs(self, file_path: Path, func_name: str, start_line: int) -> Tuple[str, int, int]:
        """Extract function definition including preceding comments."""
        with open(file_path, 'r') as f:
            lines = f.readlines()

        # Find the start of documentation (look backwards for comment block)
        doc_start = start_line - 1
        while doc_start > 0:
            line = lines[doc_start - 1].strip()
            if line.startswith('/*') or line.startswith('//') or line.startswith('*') or line == '':
                doc_start -= 1
            else:
                break

        # Find end of function (match braces)
        func_end = start_line - 1
        brace_count = 0
        in_function = False

        for i in range(start_line - 1, len(lines)):
            line = lines[i]

            # Count braces
            for char in line:
                if char == '{':
                    brace_count += 1
                    in_function = True
                elif char == '}':
                    brace_count -= 1

            # If we've closed all braces, we're done
            if in_function and brace_count == 0:
                func_end = i + 1
                break

            # If line ends with semicolon (forward declaration), that's the end
            if ';' in line and '{' not in line:
                func_end = i + 1
                break

        # Extract the function text
        function_text = ''.join(lines[doc_start:func_end])

        return function_text, doc_start + 1, func_end

    def extract_all_dead_functions(self) -> Dict[str, List[Tuple[str, str, int, int]]]:
        """Extract all dead functions with their text."""
        extracted = defaultdict(list)

        for file_path_str, funcs in self.dead_functions.items():
            file_path = Path(file_path_str)

            for func_name, line_num in funcs:
                text, start, end = self.extract_function_with_docs(file_path, func_name, line_num)
                extracted[file_path_str].append((func_name, text, start, end))

        return extracted

    def create_legacy_file(self, extracted_funcs: Dict[str, List[Tuple[str, str, int, int]]]) -> str:
        """Create a legacy file with all dead functions."""
        legacy_content = """/**
 * LEGACY DEAD FUNCTIONS
 *
 * This file contains functions that were defined but never used in the current
 * SSLS models. They have been moved here for historical reference.
 *
 * Date: 2026-02-10
 * Reason: Dead code cleanup - functions never called
 */

"""

        for file_path, funcs in sorted(extracted_funcs.items()):
            legacy_content += f"\n{'='*80}\n"
            legacy_content += f"// From: {file_path}\n"
            legacy_content += f"{'='*80}\n\n"

            for func_name, text, start, end in funcs:
                legacy_content += f"// Lines {start}-{end}: {func_name}\n"
                legacy_content += text
                legacy_content += "\n\n"

        return legacy_content

    def remove_functions_from_files(self, extracted_funcs: Dict[str, List[Tuple[str, str, int, int]]]):
        """Remove dead functions from original files."""

        for file_path_str, funcs in extracted_funcs.items():
            file_path = Path(file_path_str)

            with open(file_path, 'r') as f:
                lines = f.readlines()

            # Sort by start line in reverse order (so we can delete from bottom up)
            funcs_sorted = sorted(funcs, key=lambda x: x[2], reverse=True)

            # Remove each function
            for func_name, text, start, end in funcs_sorted:
                # Delete lines start-1 to end-1 (convert to 0-indexed)
                del lines[start-1:end]

            # Write back
            with open(file_path, 'w') as f:
                f.writelines(lines)

            print(f"✓ Removed {len(funcs)} functions from {file_path}")

    def print_summary(self, extracted_funcs: Dict[str, List[Tuple[str, str, int, int]]]):
        """Print summary of what will be moved."""
        print(f"\n{'='*80}")
        print("DEAD FUNCTION EXTRACTION SUMMARY")
        print(f"{'='*80}\n")

        total = 0
        for file_path, funcs in sorted(extracted_funcs.items()):
            print(f"{file_path}:")
            for func_name, _, start, end in funcs:
                lines_count = end - start + 1
                print(f"  - {func_name}() [lines {start}-{end}, {lines_count} lines]")
                total += lines_count
            print()

        print(f"Total: {sum(len(f) for f in extracted_funcs.values())} functions, ~{total} lines\n")


def main():
    stan_dir = Path("/mnt/code/stan")
    extractor = StanFunctionExtractor(stan_dir)

    # Extract all dead functions
    extracted = extractor.extract_all_dead_functions()

    # Print summary
    extractor.print_summary(extracted)

    # Proceed with extraction
    print("Proceeding with extraction...")
    print("  1. Creating stan/legacy/dead_functions.stan with all dead functions")
    print("  2. Removing dead functions from original files")
    print()

    # Create legacy file
    legacy_content = extractor.create_legacy_file(extracted)
    legacy_dir = stan_dir / "legacy"
    legacy_dir.mkdir(exist_ok=True)

    legacy_file = legacy_dir / "dead_functions.stan"
    with open(legacy_file, 'w') as f:
        f.write(legacy_content)

    print(f"✓ Created {legacy_file}")

    # Remove from original files
    print()
    extractor.remove_functions_from_files(extracted)

    print("\n✅ Dead functions moved to legacy!")
    print("\nNext: Verify models compile")


if __name__ == "__main__":
    main()
