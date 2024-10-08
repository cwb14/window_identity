import sys

def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    cleaned_lines = []
    for line in lines:
        stripped_line = line.strip()
        if stripped_line == '###':
            # Only add "###" if it's not already the last line added
            if cleaned_lines and cleaned_lines[-1] != '###':
                cleaned_lines.append('###')
        elif stripped_line:
            parts = stripped_line.split()
            if len(parts) >= 2:
                cleaned_line = parts[0] + '\t' + parts[1]
                cleaned_lines.append(cleaned_line)
            else:
                # Handle lines with fewer than two parts gracefully
                print(f"Warning: Line '{stripped_line}' in '{file_path}' does not have at least two parts. Skipping.")
    
    # Ensure the list starts and ends with "###"
    if not cleaned_lines:
        return ['###']
    
    if cleaned_lines[0] != '###':
        cleaned_lines.insert(0, '###')
    
    if cleaned_lines[-1] != '###':
        cleaned_lines.append('###')
    
    return cleaned_lines

def main():
    if len(sys.argv) < 2:
        print("Usage: python anchor_builder_cleaned.py <file1> [file2] [file3] ...")
        sys.exit(1)

    all_output_lines = ['###']  # Start with "###" as per requirement
    for file_path in sys.argv[1:]:
        output_lines = process_file(file_path)
        # Avoid consecutive "###" when merging outputs from multiple files
        if all_output_lines[-1] == '###' and output_lines[0] == '###':
            output_lines = output_lines[1:]
        all_output_lines += output_lines

    # Ensure the final output ends with "###"
    if all_output_lines[-1] != '###':
        all_output_lines.append('###')

    # Print the final concatenated output
    for line in all_output_lines:
        print(line)

if __name__ == "__main__":
    main()
