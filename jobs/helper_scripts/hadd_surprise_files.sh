files_ana_list=/exp/uboone/data/users/jbateman/book/Trident_EPEM/v10_04_07_20/Trident_Overlay/SURPRISE_reco2/filesana.list
output_path=/exp/uboone/data/users/jbateman/workdir/DarkNews/Trident/SURPRISE 
out_file=trident_SURPRISE_100evts.root

hist_template=reco_stage_2_hist_

# Clear the temp file list
>temp_filesana.list

# Loop through each line in the files_ana_list
while IFS= read -r line; do
    # Extract the filename from the line
    filename=$(echo "$line" | awk '{print $NF}')
    # Check that the filename contains the hist_template
    if [[ "$filename" == *"$hist_template"* ]]; then
        # Append the filename to the temp_filesana.list
        echo "$filename" >> temp_filesana.list
    fi
done < "$files_ana_list"

hadd -f $output_path/$out_file @temp_filesana.list