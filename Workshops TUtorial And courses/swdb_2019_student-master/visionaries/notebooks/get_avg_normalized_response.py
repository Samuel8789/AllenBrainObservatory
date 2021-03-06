def get_avg_normalized_response(boc, session_id, cell_specimen_id, temporal_frequency=2.0):
    ''' generate normalized average response for each grating orientation
    '''
    data_set = boc.get_ophys_experiment_data(session_id)
    
    timestamps, dff = data_set.get_dff_traces(cell_specimen_ids=[cell_specimen_id])
    dff_trace = dff[0,:]
    
    stim_table = data_set.get_stimulus_table('drifting_gratings')
    
    #Calculate the mean DF/F for each grating presentation in this stimulus
    rows = list()
    for i in range(len(stim_table)):
        new_row = {
            'orientation': stim_table.orientation[i],
            'temporal_frequency': stim_table.temporal_frequency[i],
            'max_dff': dff_trace[stim_table.start[i]:stim_table.end[i]].max()
        }
        rows.append(new_row)

    cell_response = pd.DataFrame.from_dict(rows)
    tf2_response = cell_response[cell_response.temporal_frequency==temporal_frequency]
    
    mean_dff_ori = tf2_response.groupby('orientation').mean()
    
    max_response = tf2_response['max_dff'].max()
    
    return mean_dff_ori['max_dff']/max_response

