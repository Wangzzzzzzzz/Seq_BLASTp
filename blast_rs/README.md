# README of the JSON output of the BLAST Results

## JSON File names

* Left Data Set -> query dataset

* Right Data Set -> subject dataset

* _Basically, the left dataset is compared against right dataset_

## Explanation of the descriptor:

| Item Name | Description |
| ----- | ------ |
| "score" | blast score, the higer the score, the better the match is |
| "E_Value" |  the Expected Value of hits by chance, the lower E Value is, the better the hit |
| "Query Cover" | length of aligned pair over over all length of the __query sequence__. This multiplied by "Perc_Indentity" should gives similarity of the subject w.r.t query sequence. However, this derived field is not give here, as it is not an offical result obtianed from blast. |
| "Perc_Identity" | ~~percentages of identical matches in the Aligned Sequence~~ number of identical matches over number of total aligned pairs  |
| "Num_of_Identity" | number of identical matches in the aligned portion |
| "Aligned_length" | length of the aligned portion |
| "Align_Seq_Length" | length of the subject/database sequence that result in the hit |
| "Query_Seq_Length" | length of the query sequence |
| "query_start_pos" | position of query sequence when the alignment begins |
| "aligned_query_seq" | the sequence of aligned portion of query chain |
| "aligned_match_pts" | showing the matches in the query and subject |
| "aligned_sbjct_seq" | the sequence of aligned portion of subject chain |
| "sbjct_start_pos" | position of subject sequence when the alignment begins |