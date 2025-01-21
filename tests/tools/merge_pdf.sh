#!/bin/bash

# Define keywords to search for PDF files
keywords=(
    "plain"
    "noisy_1_no_rotation"
    "noisy_2_no_rotation"
    "noisy_3_no_rotation"
    "noisy_4_no_rotation"
    "linearly_transformed"
    "rotation_noisy_1"
    "rotation_noisy_2"
    "rotation_noisy_3"
    "rotation_noisy_4"
    "permuted"
    "permuted_noisy_1"
    "permuted_noisy_2"
    "permuted_noisy_3"
    "permuted_noisy_4"
    "quantized"
    "perturbed_x0"
    "random_nan"
    "truncated"
)

output_file="merged.pdf"
declare -A keyword_to_files  # Associative array to map keywords to files
declare -A seen_files  # Associative array to track seen files (for deduplication)

# Print all PDF files for debugging
echo "Found these PDF files:"
find . -maxdepth 1 -name "summary*.pdf" -type f | while read -r file; do
    echo "  $file"
done

# Function to add files to array in natural sort order
add_files_sorted() {
    local key=$1
    local files=("${@:2}")
    if [ ${#files[@]} -gt 0 ]; then
        # Convert array to newline-separated string, sort, and store
        sorted_files=$(printf "%s\n" "${files[@]}" | sort -V)
        keyword_to_files[$key]="$sorted_files"
    fi
}

# Search for PDF files with keywords and sort them
for keyword in "${keywords[@]}"; do
    echo "Searching for keyword: $keyword"
    case "$keyword" in
        "perturbed_x0")
            files=($(find . -maxdepth 1 -type f -name "summary*perturbed_x0*.pdf"))
            add_files_sorted "$keyword" "${files[@]}"
            ;;
        "random_nan")
            files=($(find . -maxdepth 1 -type f -name "summary*random_nan*.pdf"))
            add_files_sorted "$keyword" "${files[@]}"
            ;;
        "truncated")
            files=($(find . -maxdepth 1 -type f -name "summary*truncated*.pdf"))
            add_files_sorted "$keyword" "${files[@]}"
            ;;
        *)
            files=($(find . -maxdepth 1 -type f -name "summary*${keyword}.pdf" -o -name "summary*${keyword}_*.pdf"))
            add_files_sorted "$keyword" "${files[@]}"
            ;;
    esac
done

# Clear array to store PDF files in order of keywords
pdf_files=()

# Add PDF files to the array in order of keywords, ensuring no duplicates
for keyword in "${keywords[@]}"; do
    if [[ -n "${keyword_to_files[$keyword]}" ]]; then
        while IFS= read -r file; do
            if [[ -f "$file" && -z "${seen_files[$file]}" ]]; then
                pdf_files+=("$file")
                seen_files["$file"]=1  # Mark file as seen
            fi
        done <<< "${keyword_to_files[$keyword]}"
    fi
done

# Print the array content for debugging
echo -e "\nFiles in order of keywords (no duplicates):"
printf '%s\n' "${pdf_files[@]}"

# Print total number of files found
echo -e "\nTotal files found: ${#pdf_files[@]}"

# Merge PDF files
if [[ ${#pdf_files[@]} -gt 0 ]]; then
    pdfunite "${pdf_files[@]}" "$output_file"
    echo "Merge successfully: $output_file"
else
    echo "There are no PDF files to merge."
    echo -e "\nAll PDF files in current directory:"
    find . -maxdepth 1 -name "summary*.pdf" -type f
fi