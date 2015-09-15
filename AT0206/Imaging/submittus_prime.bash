
# Submit all of the submission files!

for (( i = 10; i < 205; i++ )); do
    canfar_submit /home/ekoch/code_repos/canfar_scripts/img_pipe/archival_data/single_channel_subs/single_channel_${i}.sub ewk_casa_4_20 126e8ef0-b816-43ed-bd5f-b1d4e16fdda0
done

# for ((i = 90; i < 92; i++)) do
#     canfar_submit /home/ekoch/code_repos/canfar_scripts/img_pipe/archival_data/single_channel_subs/single_channel_${i}.sub ewk_casa_6_19 126e8ef0-b816-43ed-bd5f-b1d4e16fdda0
# done

# left=($(seq 93 96))

# running=([90]=90)

# for (( i = 91; i < 92; i++ )); do
#     running+=([$i]=$i)
# done

# ii=3

# posn=0

# while [[ $ii -ge 0 ]]; do
#     chan_ct=0
#     while [[ ${#running[@]} -ge 0 ]]; do
#         chan=${running[$chan_ct]}
#         echo "Is "${chan}" done?"
#         # Reset counter
#         if [[ ! ${running[$chan_ct]+abc} ]]; then
#             echo "Resesting channel counter."
#             chan_ct=0
#             chan=${running[$chan_ct]}
#         fi

#         if [[ -s output_files/output/single_channel_clean_${chan}.out ]]; then
#             :
#         else
#             next_chan=${left[$posn]}
#             echo "Yes. Now running "${next_chan}
#             sleep 5
#             #canfar_submit /home/ekoch/code_repos/canfar_scripts/img_pipe/archival_data/single_channel_subs/single_channel_${next_chan}.sub ewk_casa_6_19 126e8ef0-b816-43ed-bd5f-b1d4e16fdda0
#             ii=$(($ii - 1))
#             chan_ct=$(($chan_ct + 1))
#             posn=$(($posn + 1))

#             echo ${running[$chan]}
#             running[$chan]=''
#             echo ${running[$chan]}

#             running[$next_chan]=${next_chan}

#             running=( ${running[@]} )
#             echo "These channels remain running: "${running[@]}
#             echo ${#running[@]}
#         fi
#     done
# done