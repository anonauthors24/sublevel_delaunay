import os
import csv
import sys
import re

# Retrieve command-line arguments
# Directory containing the files
directory = sys.argv[1]

# Output CSV file
output_file = sys.argv[2]

# Initialize the CSV data list with headers
csv_data = [['File Name', 'Data Set', 'Function', 'Dimension', 'NumPoints', 'Complex Size', 'Delaunay Complex Size',
             'Ratio', 'Homology', 'Overall Timer', 'Initial Timer', 'Complex Timer', 'Face Timer', 'Meb Timer',
             'Graded Matrices Timer', 'Multi Chunk Timer', 'Mpfree Timer',
             'Delaunay Timer', 'Memory', 'Simplices per second', 'Time per simplex (ms)',
             'Num of mebs computed', 'Compression Rate', 'N before', 'N after',
             'Multi-chunk-Ratio']]
# Timer value pattern
timer_pattern = r'(\d+(?:\.\d+)?)'

# Iterate over the files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.txt'):
        file_path = os.path.join(directory, filename)

        # Extract values from the filename
        dataset, function, *_ = filename.split('_')

        # Read the file and extract relevant data
        with open(file_path, 'r') as file:
            file_data = file.read()

        # Parse the file data and extract the required values
        # Modify this section based on the structure of your file data

        # Extract Dimension
        dimension = None
        if 'Dimension is' in file_data:
            dimension = int(file_data.split('Dimension is ')[1].split('\n')[0])

        # Extract Points
        points = None
        if 'Read ' in file_data and ' points' in file_data:
            points = int(file_data.split('Read ')[1].split(' points')[0])

        # Extract Complex Size
        complex_size = None
        if 'Simplex tree has' in file_data:
            complex_size = int(file_data.split('Simplex tree has ')[1].split(' vertices and ')[1].split(' simplices')[0])
        
        # Extract Delaunay Complex Size
        delaunay_size = None
        if 'Total size of Del: ' in file_data:
            delaunay_size = int(file_data.split('Total size of Del: ')[1].split('\n')[0])

        # Extract Ratio
        ratio = None
        if 'Ratio: ' in file_data:
            ratio = file_data.split('Ratio: ')[1].split('\n')[0]

        # Extract Homology Dimension (if available)
        homology_dim = None
        if 'Homology dimension is' in file_data:
            homology_dim = file_data.split('Homology dimension is ')[1].split('\n')[0]
            # Convert to integer
            homology_dim = int(homology_dim)
            # Map homology dimension to corresponding label
            homology_dim = f'H_{homology_dim}'

        # Extract Timer Values
        overall_timer = None
        if 'Overall timer: ' in file_data:
            overall_timer = file_data.split('Overall timer: ')[1].split('\n')[0].strip()

        # Extract Timer Values
        overall_timer = None
        if 'Overall timer: ' in file_data:
            overall_timer = file_data.split('Overall timer: ')[1].split('\n')[0].strip()

        initial_timer = None
        if 'Inital timer: ' in file_data:
            initial_timer = file_data.split('Inital timer: ')[1].split('     (')[0].strip()
        
        complex_timer = None
        if 'Complex timer: ' in file_data:
            complex_timer = file_data.split('Complex timer: ')[1].split('     (')[0].strip()
        
        face_timer = None
        if 'Face timer: ' in file_data:
            face_timer = file_data.split('Face timer: ')[1].split('     (')[0].strip()
       
        meb_timer = None
        graded_timer = None
        chunk_timer = None
        mpfree_timer = None
        delaunay_timer = None
        
        timer_lines = file_data.split('\n')
        for line in timer_lines:
            if 'Meb timer: ' in line:
                meb_timer_value = line.split('Meb timer: ')[1].split('     (')[0].strip()
                if meb_timer_value:
                    meb_timer = float(meb_timer_value)
            elif 'Graded matrices timer: ' in line:
                graded_timer_value = line.split('Graded matrices timer: ')[1].split('     (')[0].strip()
                if graded_timer_value:
                    graded_timer = float(graded_timer_value)
            elif 'Multi chunk timer: ' in line:
                chunk_timer_value = line.split('Multi chunk timer: ')[1].strip().split()[0]
                if chunk_timer_value:
                    chunk_timer = float(chunk_timer_value)
            elif 'Mpfree timer: ' in line:
                mpfree_timer_value = line.split('Mpfree timer: ')[1].strip().split()[0]
                if mpfree_timer_value:
                    mpfree_timer = float(mpfree_timer_value)
            elif 'Delaunay timer: ' in line:
                delaunay_timer_value = line.split('Delaunay timer: ')[1].split('     (')[0].strip()
                if delaunay_timer_value:
                    delaunay_timer = float(delaunay_timer_value)        
        # Remove leading and trailing spaces from timer values
        #overall_timer = overall_timer.strip() if overall_timer is not None else None
        #initial_timer = initial_timer.strip() if initial_timer is not None else None
        #complex_timer = complex_timer.strip() if complex_timer is not None else None
        #face_timer = face_timer.strip() if face_timer is not None else None
        #meb_timer = meb_timer.strip() if meb_timer is not None else None
        #graded_timer = graded_timer.strip() if graded_timer is not None else None
        #chunk_timer = chunk_timer.strip() if chunk_timer is not None else None
        #mpfree_timer = mpfree_timer.strip() if mpfree_timer is not None else None
        #delaunay_timer = delaunay_timer.strip() if delaunay_timer is not None else None
        
        # Extract Memory
        memory = None
        if 'Memory in the end: ' in file_data:
            memory = file_data.split('Memory in the end: ')[1].split('\n')[0]

        # Extract Simplices per second
        simplices_per_second = None
        if 'Simplices per second: ' in file_data:
            simplices_per_second = file_data.split('Simplices per second: ')[1].split('\n')[0]

        # Extract Time per simplex
        time_per_simplex = None
        if 'Time per simplex (in microseconds): ' in file_data:
            time_per_simplex = file_data.split('Time per simplex (in microseconds): ')[1].split('\n')[0]


        # Extract Num of mebs computed
        num_mebs_computed = None
        if 'Computed ' in file_data and ' mebs' in file_data:
            num_mebs_computed = int(file_data.split('Computed ')[1].split(' mebs')[0])

        # Extract N before
        n_before = int(file_data.split('N before=')[1].split('\n')[0].strip()) if 'N before=' in file_data else None

        # Extract N after
        n_after = int(file_data.split('N after =')[1].split('\n')[0].strip()) if 'N after =' in file_data else None

        # Extract Compression Rate (if available)
        compression_rate = None
        if 'Compression rate: ' in file_data:
            compression_rate = file_data.split('Compression rate: ')[1].split('\n')[0]

        # Extract Multi-chunk-Ratio
        multi_chunk_ratio = None
        if 'Multi-chunk-Ratio: ' in file_data:
            multi_chunk_ratio = file_data.split('Multi-chunk-Ratio: ')[1].split('\n')[0]

        # Add the extracted values to the CSV data list
        csv_data.append([filename, dataset, function, dimension, points,
                         complex_size, delaunay_size, ratio, homology_dim,
                         overall_timer, initial_timer, complex_timer,
                         face_timer, meb_timer, graded_timer, chunk_timer,
                         mpfree_timer, delaunay_timer, memory,
                         simplices_per_second, time_per_simplex, num_mebs_computed,
                         compression_rate, n_before, n_after, multi_chunk_ratio])
# Write the CSV data to the output file
with open(output_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(csv_data)

print(f'Summaries generated and saved to {output_file}.')

