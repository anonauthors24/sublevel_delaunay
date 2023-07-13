if [ $# -eq 0 ]; then
    echo "./bifiltered_graph results_dir size"
    exit 1
fi



# Check if compute_bifiltered_graph exists
if [ ! -f compute_bifiltered_graph ]; then
    # Find C++ compiler
    if command -v g++ >/dev/null 2>&1; then
        compiler="g++"
    elif command -v clang++ >/dev/null 2>&1; then
        compiler="clang++"
    else
        echo "No C++ compiler found. Please install g++ or clang++."
        exit 1
    fi

    # Compile compute_functions.cpp
    echo "Compiling compute_bifiltered_graph.cpp .."
    $compiler --verbose -std=c++17 compute_bifiltered_graph.cpp -o  compute_bifiltered_graph 
    echo "compute_bifiltered_graph completed."
fi


TMPOUTPUT=$(mktemp /tmp/reduced_edges.XXXXX)
echo "Writing to tmp file ${TMPOUTPUT}"
#TMPFILE=results
TMPFILE=$1
size=$2
mkdir -p "$TMPFILE"/examples
mkdir -p "$TMPFILE"/bifiltered_graph
mkdir -p "$TMPFILE"/bifiltered_graph_outputs

dir_example="$TMPFILE"/examples
dir_output="$TMPFILE"/bifiltered_graph


echo "Computing examples of $size vertices" # and 2000 vertices will blow up memory
directory=./"$dir_example"

# Check if the directory exists
if [ -d "$directory" ]; then
    # Declare an empty array
    files=()

    # Iterate over each file in the directory
    for file in "$directory"/*density_${size}; do
        # Check if the current item is a file (not a directory)
        if [ -f "$file" ]; then
            # Add the file name to the array
            files+=("$(basename "$file")")

            # Check if the file name ends with "1000"
            #if [[ $filename == *_1000 ]]; then
            #    # Add the file name to the array
            #    files+=("$filename")
            #fi
        fi
    done


    # Loop through the array of file names and compute the complex
    for filename in "${files[@]}"; do
        echo "Processing... $filename"
        #echo "./"$dir_example"/${filename}"
        echo "Compute bifiltered complete graphs from ${filename}"
        ./compute_bifiltered_graph ./"$dir_example"/${filename} > ./"$dir_output"/${filename}
        ./reduce_edges ./"$dir_output"/${filename} "$TMPOUTPUT" > "$TMPFILE"/bifiltered_graph_outputs/${filename}.txt 
#        if [ "$filename" == "S2_density_2000" ]; then
#            # S2_density_2000.txt ignored
#            echo "Skipping $filename due to resource constraints"
#            continue
#        fi
        ./filtration_domination "$TMPOUTPUT" "$dir_example"/${filename} >> "$TMPFILE"/bifiltered_graph_outputs/${filename}.txt 
    done
else
    echo "Directory not found: $directory"
fi

echo "Removing ${TMPOUTPUT}"
rm "$TMPOUTPUT"
#python3 ./generate_summaries_bifiltered_graph.py ./"$TMPFILE"/bifiltered_graph_outputs/ ./"$TMPFILE"/graph_summaries.csv
#python3 ./generate_summaries.py ./"$dir_output" ./"$TMPFILE"/summaries.csv

