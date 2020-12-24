import pickle


def fill_dict(data_path, pickle_path):
    data_file = open(data_path)
    d = dict()
    #go over the data file
    #skip the first line
    next(data_file)
    for line in data_file:
        line_as_array = line.split(';')
        start, stop, chro = int(line_as_array[4][1:-1]), int(line_as_array[4][1:-1]), line_as_array[8][3:-1]
        if chro == "loroplast":
            chro = -1
        else: chro = int(chro)


        ###########################################
        ####shuold add the oposit strand as well###
        ###########################################
        d[str(line_as_array[0][1:-1])] = (start, stop, chro)
    pickle.dump(d, open(pickle_path, "wb"))
    return d

def offtarget_ininput_genes(genes_lst, genes_pos_dict, start, stop, chro):
    '''
    :param genes_lst: input of CRISPys
    :param genes_pos_dict:
    :param start: of target
    :param stop: of target
    :param chr: chromosome number of target
    :return:
    '''
    for gene in genes_lst:
        gene_cord = genes_pos_dict[gene]
        if is_overlap(gene_cord, start, stop, chro):
            return True
    return False


def is_overlap(gene_cord, start, stop, chro):
    if chro != gene_cord[2]:
        return False
    return not (start > gene_cord[1] or stop < gene_cord [0]) #chack if this is defenition true




def p_d(d):
    for item in d.items():
        print(item)







if __name__ == "__main__":
    p_d(fill_dict("D:\\LabNotCode\\Cdata\\annotation.sly.csv", "D:\\Lab\\CrispysV1.6\\genes_locations_dict.p"))