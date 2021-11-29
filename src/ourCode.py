# -*- coding: utf-8 -*-
"""
Simulation of a genome evolving inside an environment without mutations on
coding sequences.
"""
import numpy as np
import os
import re
import matplotlib.pyplot as plt
import pandas as pd

#=======================================================================
#                   READ INPUT FILES
#=======================================================================
def parse_namefile_ini(line):
    """
    Parse a line from a .ini file to get the described file location
    
    Parameters
    ----------
        line : str
            Line from the .ini file, containing the location of a file
            of interest.
    
    Returns
    -------
        processed_line : str
            Address of the file of interest
            
    >>> parse_namefile_ini("file = file/path.txt")
    'file/path.txt'
    """
    
    # get the address of the file of interest
    processed_line = line.split(" ")[2]
    # removes the newline character at the end of the string
    processed_line = processed_line.rstrip()
    return processed_line


def pos_out_genes(file_ini, folder):
    """
    Return a description of the genome described by a .ini file.
    
    Parameters
    ----------
    file_ini : str
        Initialization .ini file name.
    folder : str
        folder of the .ini file
        
    Returns
    -------
    start : list of ints
        Positions of gene starts
    end : list of ints
        Positions of gene ends
    barr : list of ints
        Positionds of topological barriers
    out : Numpy array
        2-D array of ints. Each line represents an open interval containing
        no gene nor barrier
    genome_size : int
        Number of base pairs in the genome.
    """
    
    f = open(file_ini, 'r')
    f_lines = f.readlines() # List of the file_ini lines
    file_gff = folder + parse_namefile_ini(f_lines[1]) # GFFF file address
    file_tss = folder + parse_namefile_ini(f_lines[2]) # TSS file address
    file_tts = folder + parse_namefile_ini(f_lines[3])# TTS file address
    file_barr = folder + parse_namefile_ini(f_lines[4]) # prot.dat file address
    
    ## Genome size
    fgff = open(file_gff, 'r')
    fgff_content = fgff.read() # List of the file_gff lines
    Ngen = re.findall("##sequence-region .* 1 ([0-9]+)",fgff_content)[0]
                      
    ## List of transcription start positions
    ftss = open(file_tss, 'r')
    ftss_lines = ftss.readlines() # List of the file_tss lines
    start = []
    for x in ftss_lines[1:]:
        start.append(int(x.split("\t")[2]))
    
    ## List of transcription end positions
    ftts = open(file_tts, 'r')
    ftts_lines = ftts.readlines() # List of the file_tts lines
    end = []
    for x in ftts_lines[1:]:
        end.append(int(x.split("\t")[2]))
    
    ## List of topological barriers
    fbarr = open(file_barr, 'r')
    fbarr_lines = fbarr.readlines() # List of the file_barr lines
    barr = []
    for x in fbarr_lines[1:]:
        barr.append(int(x.split("\t")[1]))
    
    ## Close files
    f.close()
    fgff.close()
    ftss.close()
    ftts.close()
    fbarr.close()
    
    ### List of open position intervals in the genome where there is 
    ### no gene
    out = pos_out_from_pos_lists(start, end, barr)
    return start, end, barr, out, int(Ngen)


def target_expression(environment_file):
    """
    Get the target expression from given file.
    
    Parameters
    ----------
    environment_file : str
        Path and name of the environment file (giving the target relative
        expression levels)

    Returns
    -------
    target_expresstion : Numpy array
        Array of floats representing target relative expression level for each
        gene.
    """
    
    with open(environment_file, "r") as environment:
        target_exp = re.findall("[0-9]+\s+([0-9\.]+)", environment.read())
    return np.asarray(target_exp, dtype=float)


#=======================================================================
#                   WRITE OUTPUT FILES
#=======================================================================
def update_files(genome_size, genes_start_pos, genes_end_pos, barriers_pos,
                 gff_file, tss_file, tts_file, barr_file):
    """Write the initialization files for the transcription simulation.
    
    The event can either be a genome inversion (with proba inversion_proba),
    or an indel.
    
    Parameters
    ----------
    genome_size : int
        Genome size in base pair.
    genes_start_pos : Numpy array
        Array of ints representing the begining position of genes.
    genes_end_pos : Numpy array
        Array of ints representing the ending position of genes.
    barriers_pos : Numpy array
        Array of ints representing the position of barriers.
    gff_file : str
        Name of the .gff file (gene positions).
    tss_file : str array
        Name of the TSS file (gene start positions).
    tts_file : str
        Name of the TTS file (gene end positions).
    barr_file : str
        Name of the .dat file containing barrier positions
        
    Note
    ----
    Nothing is returned, but the files are created or updated.
    """
    
    sequence_name = gff_file[:-4]
    new_gff = open(gff_file, "w")
    new_tss = open(tss_file, "w")
    new_tts = open(tts_file, "w")
    new_barr = open(barr_file, "w")
    ### Headers
    new_gff.writelines(["##gff-version 3\n",
                        "#!gff-spec-version 1.20\n",
                        "#!processor NCBI annotwriter\n",
                        "##sequence-region " + sequence_name + " 1 "
                        + str(genome_size) + "\n",
                        sequence_name + "\tRefSeq\tregion\t1\t"
                        + str(genome_size) + "\t.\t+\t.\tID=id0;Name="
                        + sequence_name + "\n"])
    new_tss.write("TUindex\tTUorient\tTSS_pos\tTSS_strength\n")
    new_tts.write("TUindex\tTUorient\tTTS_pos\tTTS_proba_off\n")
    new_barr.write("prot_name\tprot_pos\n")
    ### Body
    for gene_index in range(len(genes_start_pos)):
        start = genes_start_pos[gene_index]
        end = genes_end_pos[gene_index]
        size = end - start
        if size > (genome_size / 2):
            # Gene oriented "-" but crossing the origin
            orient = "-"
        elif size > 0:
            orient = "+"
        elif size > (-genome_size / 2):
            orient = "-"
        else:
            # Gene oriented "+" but crossing the origin
            orient = "+"
        new_gff.write(sequence_name + "\tRefSeq\tgene\t" + str(start) + "\t"
                      + str(end) + "\t.\t" + orient + "\t.\tID=g1;Name=g"
                      + str(gene_index + 1) + "\n")
        new_tss.write(str(gene_index) + "\t" + orient + "\t" + str(start)
                      + "\t.2\n")
        new_tts.write(str(gene_index) + "\t" + orient + "\t" + str(end)
                      + "\t1.\n")
    for barrier in barriers_pos:
        new_barr.write("hns\t" +str(barrier) + "\n")
    ### Close
    for file in [new_gff, new_tss, new_tts, new_barr]:
        file.close()


def copy_genome(PARAMS):
    """
    Copy the last accepted genome to save it
    
    Parameters
    ----------
    PARAMS : Python list
        list of the following parameters:
            NEXT_GEN_GFF : str
                Name of the .gff file (gene positions) at the next generation.
            NEXT_GEN_TSS : str
                Name of the TSS file (gene start positions) 
                at the next generation.
            NEXT_GEN_TTS : str
                Name of the TTS file (gene end positions)
                at the next generation.
            NEXT_GEN_BARRIERS : str
                Name of the .dat file containing barrier positions 
                at the next generation.
            LAST_ACCEPTED_GENOME : list of str
                Name of the equivalent files, where we want to copy the previous ones.
    Note
    ----
    Nothing is returned, but the files are created or updated.
    """
    
    for k, filename in enumerate(PARAMS[:-1]):
        with open(filename, "r") as read_file:
            with open(PARAMS[-1][k], "w") as write_file:
                write_file.write(read_file.read())
  
  
#=======================================================================
#                   PROCESS THE GENOME
#=======================================================================
def pos_out_from_pos_lists(genes_start_pos, genes_end_pos, barriers_pos):        
    """Generate the list of postition intervals outside barriers and genes.
    
    Parameters
    ----------
    genes_start_pos : Numpy array
        Array of ints representing the begining position of genes.
    genes_end_pos : Numpy array
        Array of ints representing the ending position of genes.
    barriers_pos : Numpy array
        Array of ints representing the position of barriers.
    
    Returns
    -------
    out_positions : Numpy array
        2-D array of ints. Each line represents an open interval containing
        no gene nor barrier.range
    """
    
    limits = np.sort(np.hstack((genes_start_pos, genes_end_pos, barriers_pos,
                                barriers_pos)))
    # First interval is from the last coding position to the first one
    # as the genome is circular
    out_positions = np.array([limits[-1], limits[0]])
    for k in range(1, len(limits)-1, 2):
        out_positions = np.vstack((out_positions, [limits[k], limits[k+1]]))
    return out_positions


def expression_simulation(params_file, out_file, gene_start_pos):
    """Run  the expression  simulation with given parameters.
    
    Parameters
    ----------
    params_file : str
        Path and name of the parameters file.
    out_file : str
        Path and name of the output file.
    gene_start_pos : Numpy array
        Array of ints representing the begining position of genes.
        
    Returns
    -------
    transcript_numbers : Numpy array
        Array of ints representing the number of transcripts  for each gene,
        ordered by gene ID.
        
    Note
    ----
    Also generate a text file with the given out_file name.
    """
    
    # Execute the command line to run the simulation in a terminal
    term_command = ("python3 TwisTranscripT/start_simulation.py " + params_file + " > "
                    + out_file)
    os.system(term_command)
    # Open the output file of previous command to get the expression
    # profile
    with open(out_file, "r") as out:
        file_content = out.read()
        transcript_nbs = re.findall("Transcript ID [0-9]+ : ([0-9]+)",
                                    file_content)
        transcript_nbs = np.asarray(transcript_nbs, dtype=int)
        # Transcripts are reindexed, we fetch their starting position to map
        # transcription values to the right gene.
        # (This issue was fixed in the last version of 
        # start_simulation.py by its author. This code is compatible
        # with both versions)
        starts = re.findall("\n[0-9] +[0-9] +[0-9] +[\-0-9]+ +([0-9]+)",
                            file_content)
        starts = np.asarray(starts, dtype=int)
    # Reorder transcription values
    transcription_values = np.full(len(starts), np.nan)
    for k, start_pos in enumerate(gene_start_pos):
        transcript_id, = np.where(starts == start_pos)
        transcription_values[k] = transcript_nbs[transcript_id[0]]
    return transcription_values


def compute_fitness(observed_transcript_numbers, target_frequencies):
    """Compute the fitness of an individual with given gene expression pattern
    in given environment.
    
    Parameters
    ----------
    observed_transcript_numbers : Numpy array
        Array of ints representing the observed number of transcripts for each
        gene of the invidual.
    target_frequencies : Numpy array
        Array of floats representing target relative expression level for each
        gene (environment).
    
    Returns
    -------
    fitnness : float
        computed fitness of the invidual, following the formula:
            fitness = exp(-sum(|ln(observed_for_gene_i / target_for_gene_i)|))
    """
    
    observed_frequencies = (observed_transcript_numbers
                            / np.sum(observed_transcript_numbers))
    ln_freqs = np.log(observed_frequencies / target_frequencies)
    return np.exp(-np.sum(np.abs(ln_freqs)))


#=======================================================================
#                   GENERATE MUTATIONS
#=======================================================================
def sample(out, Ngen, u):
    """
    Sample a location from a given list of intervals that satisfies the "distance to bounds condition" 
    (there is at least u nucleotides between the right and left bound and the sampled mutation position)
    This allows a deletion to happen without being to close to a gene or barrier.
    
    Parameters
    ----------
    out : Python list
        list of the open position intervals in the genome where there are not any genes   
    Ngen : int
        the length of the genome
    u : int
        unit of length of nucleotides.
        
    Returns
    -------
    mut_pos : int
        sampled location of the mutation
    """
    
    space = False # Used to make sure there is enough space between the right and left bound and the sampled mutation position
    i=100 # Max number of iterations
    
    while (not space):
        # Repeat until we find a mutation position that satisfies the distance to the bounds condition
        # or until we have repeated this loop a 100 times
        ### Sample the interval
        ## For each interval, calculate the probability to sample from it
        N = 0 # Total length of the intervals
        intlen = [] # List of the cummulative length of each intervals
        for x in out:
            a = x[0] # Left bound of the interval xout_positions[i][1]
            b = x[1] # Right boud of x
            if a>b:
                # Then x is the last interval of the plasmid with b after the first position of the genome
                xlen = abs(Ngen - a + b - 1)
                N += xlen
                intlen.append(N)
            else:
                xlen = b-a-1 # Length of the interval
                N += xlen
                intlen.append(N)
        probas = np.array(intlen)/N # List of the cummulative probabilities to sample in each interval
        ## Sample the interval
        p = np.random.uniform(0,1) # Draw a random number between 0 and 1
        sint = min(np.where(p<probas)[0]) # Index of the sampled interval
    
        ### Sample the exact location of the mutation
        a = out[sint][0]
        b = out[sint][1]

        if a>b:
            # Then the selected interval is the last interval of the plasmid with b after the first position of the genome
            mut_pos = np.random.randint(a+1, Ngen+b) # Sample the location of the mutation
            if mut_pos > Ngen:
                # Then the sampled position is located after the first position of the plasmid
                
                if (mut_pos <= Ngen + out[sint][1] - u) and (mut_pos >= out[sint][0] + u):
                    # There is enough space between the right and left bound and the sampled mutation position
                    space = True
                    mut_pos = mut_pos - Ngen
                   
        else:
            mut_pos = np.random.randint(a+1, b) # Sample the location of the mutation
            if (mut_pos <= out[sint][1] - u) and (mut_pos >= out[sint][0] + u):
                # There is enough space between the right and left bound and the sampled mutation position
                space = True
            
        i-=1
        if i==0:
            raise RuntimeError("A mutation position that that satisfies the distance to the bounds condition was not found")
    
    return mut_pos


def evolutive_event(discret_step, inversion_proba, genome_size,
                    genes_start_pos, genes_end_pos, barriers_pos,
                    out_positions, p_insertion):
    """Generate an evolutive event on given genome.
    
    The event can either be a genome inversion (with probability inversion_proba),
    or an indel. An indel is an insertion with probability p_insertion
    
    Parameters
    ----------
    inversion_proba : float
        Probability for the event to be an inversion.
    genome_size : int
        Genome size in base pair.
    genes_start_pos : Numpy array
        Array of ints representing the begining position of genes.
    genes_end_pos : Numpy array
        Array of ints representing the ending position of genes.
    barriers_pos : Numpy array
        Array of ints representing the position of barriers.
    out_positions : Numpy array
        2-D array of ints. Each line represents an open interval containing
        no gene nor barrier.
    p_insertion : float
        Probability for an indel event to be an insertion.
    Returns
    -------
    genome_size : ine
        Updated value of genome_size after the event.
    genes_start_pos : Numpy array
        Updated value of genes_start_pos after the event.
    genes_end_pos : Numpy array
        Updated value of genes_end_pos after the event.
    barriers_pos : Numpy array
        Updated value of barriers_pos after the event.
    genome_size: int
        Updated value of the genome size.    
    """
    
    # Draw a position for the mutation
    event_position = sample(out_positions, genome_size, discret_step)
    if np.random.rand() < inversion_proba:
        # The event will be an inversion
        # Draw the second position for the mutation
        event_position2 = sample(out_positions, genome_size, discret_step)
        return genome_inversion(genome_size, genes_start_pos, genes_end_pos,
                                barriers_pos, min(event_position,
                                                  event_position2),
                                max(event_position, event_position2))
    else:
        return indel(discret_step, genome_size, genes_start_pos, genes_end_pos,
                     barriers_pos, out_positions, p_insertion)


def indel(u, genome_size, genes_start_pos, genes_end_pos, barriers_pos, out_positions, p_insertion):
    """
    Delete or insert in the plasmid a unit with length u in base pairs 
    
    Parameters
    ----------
    u : int
        Unit of length in base pairs that is deleted or inserted.
    genome_size : int
        Genome size in base pairs.
    genes_start_pos : Numpy array
        Array of ints representing the begining position of genes.
    genes_end_pos : Numpy array
        Array of ints representing the ending position of genes.
    barriers_pos : Numpy array
        Array of ints representing the position of barriers.
    out_positions : Numpy array
        2-D array of ints. Each line represents an open interval containing
        no gene nor barrier.
    p_insertion : float
        Probability of the event to be an insertion and not a deletion.
        
    Returns
    -------
    event_type : str
        Equal to "insertion" or "deletion".
    new_genes_start_pos : Numpy array
        Updated value of genes_start_pos after inversion.
    new_genes_end_pos : Numpy array
        Updated value of genes_end_pos after inversion.
    new_barriers_pos : Numpy array
        Updated value of barriers_pos after inversion.
    genome_size: int
        Updated value of the genome size.
    """
    
    ### Sample the indel position
    indel_pos = sample(out_positions, genome_size, u)
    
    ### Initialize the new positions
    new_genes_start_pos = np.copy(genes_start_pos)
    new_genes_end_pos = np.copy(genes_end_pos)
    new_barriers_pos = np.copy(barriers_pos)
    
    ### Choose whether it is an insertion or a deletion
    p = np.random.uniform(0,1) # Draw a random number between 0 and 1
    if p<p_insertion:
        # It is an insertion
        event_type = "insertion"
        for i in range( len(genes_start_pos) ):
            if genes_start_pos[i] >= indel_pos :
                # Update gene start positions
                new_genes_start_pos[i] += u
            if genes_end_pos[i] >= indel_pos :
                # Update gene end positions
                new_genes_end_pos[i] += u
            if barriers_pos[i] >= indel_pos :
                # Update barrieres positions
                new_barriers_pos[i] += u
        # Update genome size
        genome_size += u
        
    else:
        # It is a deletion
        event_type = "deletion"        
        for i in range( len(genes_start_pos) ):        
            if genes_start_pos[i] >= indel_pos :
                # Update gene start positions
                new_genes_start_pos[i] -= u
            if genes_end_pos[i] >= indel_pos :
                # Update gene end positions
                new_genes_end_pos[i] -= u
            if barriers_pos[i] >= indel_pos :
                # Update barrieres positions
                new_barriers_pos[i] -= u
        # Update genome size
        genome_size -= u
    
    return (event_type, genome_size, new_genes_start_pos, new_genes_end_pos,
            new_barriers_pos)


def genome_inversion(genome_size, genes_start_pos, genes_end_pos, barriers_pos,
                     inversion_start, inversion_end):        
    """Perform a genome innversion on given genome at given positions.
    
    Parameters
    ----------
    genome_size : int
        Genome size in base pair.
    genes_start_pos : Numpy array
        Array of ints representing the begining position of genes.
    genes_end_pos : Numpy array
        Array of ints representing the ending position of genes.
    barriers_pos : Numpy array
        Array of ints representing the position of barriers.
    inversion_start : int
        Position of the beginning of the inversion in the genome.
    inversion_end : int
        Position of the end of the inversion in the genome.
    
    Returns
    -------
    event_type : str
        Always equal to "inversion"
    genome_size : int
        Unchanged value of genome_size.
    genes_start_pos : Numpy array
        Updated value of genes_start_pos after inversion.
    genes_end_pos : Numpy array
        Updated value of genes_end_pos after inversion.
    barriers_pos : Numpy array
        Updated value of barriers_pos after inversion.
        
    Notes
    -----
    inversion_start and inversion_end must not be inside a gene.
    inversion_end must be greater than inversion_start.
    """
    
    # Copy lists to avoid modifying them
    new_genes_start_pos = np.copy(genes_start_pos)
    new_genes_end_pos = np.copy(genes_end_pos)
    new_barriers_pos = np.copy(barriers_pos)
    # Update them
    for new_array in [new_genes_start_pos, new_genes_end_pos,
                      new_barriers_pos]:
        for (k, position) in enumerate(new_array):
            # Perform the inversion
            if (position > inversion_start) and (position < inversion_end):
                new_array[k] = inversion_start + (inversion_end - position)
    return ("inversion", genome_size, new_genes_start_pos, new_genes_end_pos,
            new_barriers_pos)


#=======================================================================
#                   SIMULATE EVOLUTION
#=======================================================================
def accept_mutation(previous_fitness, new_fitness, q):
    """
    Accept or reject the mutation, based on fitnesses comparison.
    
    A fitness increase is always accepted. A fitness loss is accepted with
    a probability exp(fitness difference / q).
    
    Parameters
    ----------
    previous_fitness : float
        Fitness of the previous generation (before mutation).
    new_fitness : float
        Fitness of the new generation (after mutation).
    q : float
        Parameters of the Monte Carlo Metropolis algorithm, controlling the
        range of accepted fitness losses.
        
    Returns
    -------
    is_accepted : bool
        True if the mutation is accepted, False else.
        
    >>> accept_mutation(0, float('Inf'), 1)
    True
    """
    
    fitness_diff = new_fitness - previous_fitness
    if fitness_diff > 0:
        return True
    else:
        return (np.random.rand() < np.exp(fitness_diff/q))


def evolution(start, end, barr, out, genome_size, initial_expression, previous_fitness, target_freqs, discret_step, q, inversion_proba, p_insertion, nb_generations, PARAMS):
    """
    Simulate the evolution with a Monte-Carlo Metropolis algorithm. 
    
    Parameters
    ----------
    start : Numpy array
        Array of ints representing the begining position of genes.
    end : Numpy array
        Array of ints representing the ending position of genes.
    barr : Numpy array
        Array of ints representing the position of barriers.
    out : Python list
        list of the open position intervals in the genome where there 
        are not any genes   
    genome_size : int
        Genome size in base pair.
    initial_expression : Numpy array
        Array of ints representing the initial number of transcripts 
        for each gene, ordered by gene ID.
    previous_fitness : float
        computed initial fitness of the invidual, following the formula:
            fitness = exp(-sum(|ln(observed_for_gene_i / target_for_gene_i)|))
    target_freqs : Numpy array
        Array of floats representing target relative expression level for each gene.
    discret_step : int
        Size of an indel event (in base pairs)
    inversion_proba : float
        Probability for an event to be an inversion.
    p_insertion : float
        Probability for an indel event to be an insertion.
    nb_generations : int
        Number of generations to run the simulation.
                
    PARAMS : Python list
        list of the following parameters:
            NEXT_GEN_PARAMS : str
                Path and name of the parameter file necessary for 
                the expression_simulation function, at the next generation.
            NEXT_GEN_GFF : str
                Name of the .gff file (gene positions) at the next generation.
            NEXT_GEN_TSS : str
                Name of the TSS file (gene start positions) 
                at the next generation.
            NEXT_GEN_TTS : str
                Name of the TTS file (gene end positions)
                at the next generation.
            NEXT_GEN_BARRIERS : str
                Name of the .dat file containing barrier positions 
                at the next generation.
            LAST_ACCEPTED_GENOME : list of str
                5 equivalent files for the last genome to save
                
    Return
    ----
    accepted_fitnesses : list of floats
        Evolution of the system fitness (only fixed mutations)
    proposed_fitnesses : list of floats
        Fitness following the proposed mutation (whether it is fixed or not)
    accepted_status : list of str
        For each generation, either "accepted" or "rejected" depending
        on the status of the proposed mutation
    all_types : list of str
        For each generation, the mutation type: "insertion", "deletion"
        or "inversion".
    generation_numbers : list of ints
        Number of each generation
    final_expression : list of ints
        Number of copies of each gene's transcript following the last 
        expression simulation of the accepted genome.
    """
    
    # Initialize results lists
    accepted_fitnesses = [previous_fitness]
    proposed_fitnesses = [previous_fitness]
    accepted_status = ["accepted"]
    all_types = ["initial"]
    
    generation_numbers = range(nb_generations+1)
    for generation in generation_numbers[1:]:
        # Random evolutive event
        event_type, new_size, new_start, new_end, new_barr = (
                evolutive_event(discret_step, inversion_proba, genome_size, start, end,
                                barr, out, p_insertion))
        # Update parameter files and run expression simulation.
        update_files(new_size, new_start, new_end, new_barr, PARAMS[1],
                     PARAMS[2], PARAMS[3], PARAMS[4])
        # Simulate expression
        error = True
        while(error):
            try:
                new_expression = expression_simulation(PARAMS[0], "out.txt", new_start)
                new_fitness = compute_fitness(new_expression, target_freqs)
                error = False
            except:
                pass
        # Accept or reject the mutation.
        print("Generation ", end="")
        print(generation, end=":\n")
        print(event_type + " event")
        print("Fitness: ", end="")
        print(new_fitness)        
        if accept_mutation(previous_fitness, new_fitness, q):
            final_expression = new_expression
            accepted_status.append("accepted")
            previous_fitness = new_fitness
            genome_size, start, end, barr = new_size, new_start, new_end, new_barr
            out = pos_out_from_pos_lists(start, end, barr)
            copy_genome(PARAMS[1:])
        else:
            accepted_status.append("rejected")
            
        # Keep track of each event
        accepted_fitnesses.append(previous_fitness)
        proposed_fitnesses.append(new_fitness)
        all_types.append(event_type)
    
    return(accepted_fitnesses, proposed_fitnesses, accepted_status, all_types,
           generation_numbers, final_expression)
              
if __name__=="__main__":
    # Input and output files
    # Initial parameter files
    PARAM_FOLDER = "TwisTranscripT/"
    INITIAL_PARAMETERS = PARAM_FOLDER + "params.ini" # name of the file containing the initial parameter values necessary for the expression_simulation function
    ENVIRONMENT = "environment.dat"
    # Intermediate generations parameter files
    NEXT_GEN_PARAMS = PARAM_FOLDER + "params_nextGen.ini"
    NEXT_GEN_GFF = PARAM_FOLDER + "nextGen/nextGen.gff"
    NEXT_GEN_TSS = PARAM_FOLDER + "nextGen/nextGenTSS.dat"
    NEXT_GEN_TTS = PARAM_FOLDER + "nextGen/nextGenTTS.dat"
    NEXT_GEN_BARRIERS = PARAM_FOLDER + "nextGen/nextGenProt.dat"
    # Last generation parameter files
    LAST_ACCEPTED_GENOME = [PARAM_FOLDER + "nextGen/last.gff",PARAM_FOLDER + "nextGen/lastTSS.dat",PARAM_FOLDER + "nextGen/lastTTS.dat",PARAM_FOLDER + "nextGen/lastProt.dat"]
    # Parameters to input in the evolution function
    PARAMS = [NEXT_GEN_PARAMS, NEXT_GEN_GFF, NEXT_GEN_TSS, NEXT_GEN_TTS, NEXT_GEN_BARRIERS, LAST_ACCEPTED_GENOME]
    # Color map for final plotting
    COLORS = {"initial" : "black", "deletion" : "red", "insertion" : "green", "inversion" : "purple"}
    
    # Inputs for the evolution function
    # Recommended values are indicated in comments
    target_freqs = target_expression(ENVIRONMENT)
    output_filename = input("path and name of the output csv file: ") # out.csv
    discret_step = int(input("Unit of length in base pairs that is deleted or inserted: ")) # 60
    inversion_proba = float(input("Probability for an evolutive event to be an inversion: ")) # 0.5
    p_insertion = float(input("Probability for an indel event to be an insertion: ")) # 0.5
    nb_generations = int(input("Number of generations: ")) # 30
    q = float(input("Value of q : ")) # 0.00002

    # Process the initial genome
    start, end, barr, out, size = pos_out_genes(INITIAL_PARAMETERS, PARAM_FOLDER)
    initial_expression = expression_simulation(INITIAL_PARAMETERS, "out.txt", start)
    previous_fitness = compute_fitness(initial_expression, target_freqs)

    # Simulate evolution
    (accepted_fitnesses, proposed_fitnesses, accepted_status, all_types,
     generation_numbers, final_expression) = evolution(start, end, barr, out, size,
                                     initial_expression, previous_fitness, target_freqs, discret_step, q, inversion_proba, p_insertion, nb_generations, PARAMS)
                                     
    # Print and plot results
    print("initial expression", initial_expression)
    print("final expression", final_expression)

    plt.ylim(.9*min(proposed_fitnesses), 1.1*max(accepted_fitnesses))
    plt.plot(accepted_fitnesses, linestyle="--", markersize=0, color="k", zorder=1)
    plt.scatter(generation_numbers, accepted_fitnesses, alpha=1,
                c=[COLORS[event_type] for event_type in all_types], zorder=2)
    plt.scatter(generation_numbers, proposed_fitnesses, marker="+",
                c=[COLORS[event_type] for event_type in all_types])
    plt.show()

    # Save the results
    if output_filename:
        output_matrix = np.vstack((accepted_fitnesses,
                                   proposed_fitnesses, accepted_status, all_types))
        output_df = pd.DataFrame(data=output_matrix, columns=generation_numbers,
                                 index=["system fitness", "proposed fitness",
                                        "accepted", "event type"])
        output_df.to_csv(output_filename)
