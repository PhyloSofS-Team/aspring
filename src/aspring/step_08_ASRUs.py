import numpy as np
import pandas as pd
import csv
import sys
import itertools
import argparse
import networkx as nx

from aspring import __version__


def parse_args(args):
    parser = get_arg_parser()
    return parser.parse_args(args)

def get_arg_parser():
    parser = argparse.ArgumentParser(
        description=
        'Identifies the Alternative Splicing Repetitive Units (ASRUs) on the gene.'
    )
    parser.add_argument('--gene',
                        dest='geneName',
                        type=str,
                        required=True,
                        help='name of gene')
    parser.add_argument('--dataPATH',
                        type=str,
                        required=True,
                        help='path to dir containing Thoraxe outputs')
    parser.add_argument('--version',
                        action='version',
                        version=f'aspring {__version__}')
    return parser


def sexSize(path_data, gene,
            sex):  #info sur le sex, sa taille est dans la colonne 3

    sexSizeEvents_path = f'{path_data}/data/{gene}/{gene}_sexSizeEvents.csv'
    df_size = pd.read_csv(sexSizeEvents_path)
    df_size_gene = df_size.loc[df_size['gene'] == gene].loc[df_size['sexon'] ==
                                                            sex].to_numpy()
    return df_size_gene


def eventsDup_df(path_data, gene):
    #fonction qui return la portion dans eventsDup du gene targeté
    path = f'{path_data}/data/{gene}/{gene}_eventsDup_withCols.txt'
    df = pd.read_csv(path)  #csv to pandas.dataframe
    df = df[df['typePair'] != 'No']
    df_gene = df.loc[df['gene'] == gene]
    return df_gene


def ases_df(path_data, gene):
    #fait le dataframe de l'ases du gene targete où on isole les colonnes
    # d'intérets (chemin can, alt, ...) voir 2eme variable

    ases_gene_path = f'{path_data}/{gene}/thoraxe/ases_table.csv'
    df_ases_gene = pd.read_csv(ases_gene_path)[[
        'CanonicalPath', 'AlternativePath', 'ASE', 'MutualExclusiveCanonical',
        'MutualExclusiveAlternative'
    ]]  #ICI
    df_ases_gene_np = df_ases_gene.to_numpy()

    return df_ases_gene_np


def sim_graph_initial(path_data, gene):
    #a pair can appear multiple times, select it only once for the creation of its edge
    tmp_G = eventsDup_df(path_data, gene).sort_values("rank").drop_duplicates(
        subset=['sexA', 'sexB'], keep='first')
    tmp_G = nx.from_edgelist(tmp_G.to_numpy()[:, 1:3], create_using=nx.Graph)
    return tmp_G


def sim_graph_workframe(path_data, gene):

    df_gene = eventsDup_df(path_data, gene)
    G = nx.Graph()
    for pair in df_gene.to_numpy():
        G.add_node(pair[1], length=pair[12])
        G.add_node(pair[2], length=pair[13])
        if G.has_edge(pair[1], pair[2]):
            G[pair[1]][pair[2]]['info'][0].append(pair[3])
        else:
            #edge attribut : [[evenements],statusA,statusB,colA,colB]
            G.add_edges_from([(pair[1], pair[2], {
                'info': [[pair[3]], pair[5], pair[6], pair[-2], pair[-1]]
            })])
    return G


def sim_graph_initial_bis(path_data, gene):
    df_gene = eventsDup_df(path_data, gene)
    G_bis = nx.Graph()
    for pair in df_gene.to_numpy():
        G_bis.add_node(pair[1])
        G_bis.add_node(pair[2])
        if G_bis.has_edge(pair[1], pair[2]):
            G_bis[pair[1]][pair[2]]['eve'].append(pair[3])
        else:
            #edge attribut : [[evenements],statusA,statusB,colA,colB]
            G_bis.add_edges_from([(pair[1], pair[2], {'eve': [pair[3]]})])
    return G_bis


def event_to_pairs(path_data, G, gene):  #G is sim_graph_workframe

    #{ (paire) : ([b_1,...,b_p],classe), ...} dico où les clés sont les paires
    #les valeurs sont les évènements dans lesquels cette paire est impliquée et la classe (MEX,ALT,REL,UNREL) de cette paire

    dico = {}

    for pair in G.edges:
        for eve in G[pair[0]][pair[1]]['info'][0]:
            if eve in dico:
                dico[eve].append(pair)
            else:
                dico[eve] = [pair]

    for eve in dico:
        dico[eve].append(ases_df(path_data, gene)[eve - 1])

    return dico


def get_alnCols(G, A, B, eveDups):
    tmp = eveDups.loc[(eveDups['sexA'] == A)
                      & (eveDups['sexB'] == B)].to_numpy()
    if len(tmp) > 0:
        tmp = tmp[0]
        colA, colB = G[A][B]['info'][-2].split('-'), G[A][B]['info'][-1].split(
            '-')
    else:
        colA, colB = G[A][B]['info'][-1].split('-'), G[A][B]['info'][-2].split(
            '-')
    return colA, colB


def extension_marge_pair(G, A, B, eveDups):  #G is sim_graph_workframe(gene)

    #calcul pour un alignement entre A et B si je peux étendre A ou B en N ou C terminal

    #A est dans le can et B dans le alt (???????)
    #FIXED:---!!!WARNING!!!--- : G[A][B]['info']  or G[B][A]['info'] return same thing yet colA , colB = tmp[-2].split('-'), tmp[-1].split('-') won't be the same !!!!!!!!
    tmp = G[A][B]['info']
    colA, colB = get_alnCols(G, A, B, eveDups)
    taille_A, taille_B = int(G.nodes[A]['length']), int(G.nodes[B]['length'])

    coverage_A = np.abs(int(colA[1]) - int(colA[0]) + 1) * 100 / taille_A
    coverage_B = np.abs(int(colB[1]) - int(colB[0]) + 1) * 100 / taille_B

    if coverage_A < 100 or coverage_B < 100:
        if int(colB[0]) == 1:
            marge_B_Nter = 0
            if int(colA[0]) == 1:
                marge_A_Nter = 0
            if int(colA[0]) > 1:
                marge_A_Nter = int(colA[0]) - 1
        if int(colB[0]) > 1:
            marge_B_Nter = int(colB[0]) - 1
            if int(colA[0]) == 1:
                marge_A_Nter = 0
            if int(colA[0]) > 1:
                marge_A_Nter = int(colA[0]) - 1
        marge_A_Cter = taille_A - int(colA[1])
        marge_B_Cter = taille_B - int(colB[1])
    else:
        marge_B_Nter, marge_B_Cter, marge_A_Nter, marge_A_Cter = 0, 0, 0, 0

        #si marge_B est positive alors on cest A quon etend sinon cest B

    return [marge_B_Nter, marge_B_Cter, marge_A_Nter, marge_A_Cter]


def codage(eve2pair, C, d_X):

    #C est le candidat pour faire l'extension
    #voir pdf du rapport pour les détails du codage

    isgood = True
    for evebis in eve2pair:  #on boucle sur les évènements pour vérifier si l'extension est compatible avec tout les évènements
        status_ext = (C in eve2pair[evebis][-1][0].split('/'), C
                      in eve2pair[evebis][-1][1].split('/'),
                      int(C in eve2pair[evebis][-1][0].split('/')) +
                      int(C in eve2pair[evebis][-1][1].split('/')))
        print('-----------------------------', C,
              eve2pair[evebis][-1][0].split('/'),
              eve2pair[evebis][-1][1].split('/'), status_ext, d_X[evebis],
              evebis)
        if d_X[evebis][2] == 2 and status_ext[2] == 1:
            isgood = False
            break
        if d_X[evebis][2] == 0 and status_ext[2] == 1:
            isgood = False
            break
        if d_X[evebis][2] == 1 and status_ext != d_X[evebis]:
            isgood = False
            break

    return isgood


def extensionbis(path_data, gene):

    allpath_path = f'{path_data}/data/{gene}/{gene}_canonical_path.txt'
    allPaths = {}  # une ligne = gene + chemin can entier
    with open(allpath_path, 'r') as f:
        for lines in f.readlines():
            temp = lines.rstrip('\n').split()
            allPaths[temp[0]] = temp[1]

    g = sim_graph_workframe(path_data, gene)

    #creation du squelette du sim graph (squelette car aucun extension n'est faite)
    eve2pair = event_to_pairs(path_data, g,
                              gene)  #dico des evenements aux paires
    eveDup = eventsDup_df(path_data, gene)

    mem_nodes = {}

    for pair in g.edges():  #je loop sur toutes les paires
        print(
            '----------------------------------------------------------------------------PAIR',
            pair)
        A, B = pair[0], pair[1]
        d_A = {}
        d_B = {}
        if A in mem_nodes:
            A = mem_nodes[A]
        if B in mem_nodes:
            B = mem_nodes[B]
        print('mem_nodes', mem_nodes)
        print('HELL-o', pair, A, B)
        #if not g.has_edge(A,B):  #WARNING : temp solution
        #	continue
        extension_marge = extension_marge_pair(
            g, A, B, eveDup)  #je calcul si je peux etendre
        print('extension_marge', extension_marge, A, B)

        if extension_marge[0] > 0 or extension_marge[1] > 0:  #Si j'étend A
            for eve in eve2pair:
                d_A[eve] = (A in eve2pair[eve][-1][0].split('/'), A
                            in eve2pair[eve][-1][1].split('/'),
                            int(A in eve2pair[eve][-1][0].split('/')) +
                            int(A in eve2pair[eve][-1][1].split('/')))
        if extension_marge[2] > 0 or extension_marge[3] > 0:  #si j'étend B
            for eve in eve2pair:
                d_B[eve] = (B in eve2pair[eve][-1][0].split('/'), B
                            in eve2pair[eve][-1][1].split('/'),
                            int(B in eve2pair[eve][-1][0].split('/')) +
                            int(B in eve2pair[eve][-1][1].split('/')))
        flag = False
        for ind in np.where(
                np.array(extension_marge) > 0
        )[0]:  #je considère l'extension pour toutes les marges non nulles
            marge = extension_marge[ind]
            print('MARGE', marge)
            for eve in g.edges[(
                    A, B
            )]['info'][0]:  #pour chaque paire on boucle sur ses évènements
                print(eve, pair)
                #attention si liste contenient plusieurs fois même valeur, index() renvoie juste indice de premiere occurence
                if ind <= 1:
                    X = A  #on etend A
                    d_X = d_A
                    if ind == 0:
                        direction = 'Nter'
                    else:
                        direction = 'Cter'
                else:
                    X = B  #on etend B
                    d_X = d_B
                    if ind == 2:
                        direction = 'Nter'
                    else:
                        direction = 'Cter'

                if ind <= 1:
                    p = 1
                else:
                    p = 2
                #WARNING!!problem here, p is info relative to event but this event info not used yet
                if g.edges[(A, B)]['info'][
                        p] == 'can':  #path_bool me dit si je regarde le chemin can ou alt
                    path_bool = 0  #chemin can
                else:
                    path_bool = 1  #chemin alt
                #X est celui qu'on étend
                if flag == False:
                    path = eve2pair[eve][-1][path_bool].split(
                        '/')  #je recupere le chemin can ou alt
                print('NOT UPDATED PATH ???', path, X)
                if X not in path:  #si X n'est pas dans le chemin can ou alt, cad que c'est une ancre
                    path = allPaths[gene][1].split(
                        '/')  #je recupere le chemin canonique global
                    if X not in path:
                        continue

                print('PATH', path)
                while marge > 0:
                    #on vérifie la taille du chemin et qu'il y a bien des éléments avant ou après X
                    #C est le candidat pour faire l'extension (qui sera de la forme C/X ou X/C)
                    print('before', marge)
                    print('PATH IN WHILE', path)
                    print('X', X)
                    if len(path) > 1 and (
                        (direction == 'Cter' and path.index(X) < len(path) - 1)
                            or (direction == 'Nter' and path.index(X) > 1)):
                        if direction == 'Cter':
                            C = path[
                                path.index(X) +
                                1]  # en Cter le candidat est le noeud pointe par X
                        else:
                            C = path[
                                path.index(X) -
                                1]  #en Nter le candidat est le noeud qui pointe X

                        if C == 'start' or C == 'stop' or (
                                C
                                not in eve2pair[eve][-1][path_bool].split('/')
                        ) or C.split('_')[0] == '0' or path[-1] == C or path[
                                0] == C:  #check if C is not start or stop node of ESG
                            break
                        else:
                            if C.split(
                                    '_'
                            )[0] == '0':  #check if C is not a node '0_blabla' which are not real sex
                                break

                            else:
                                print('I AM C', C)
                                if '$' in C:
                                    print('AAAAAAAAAAAAAAAAAAAAAAAAAAA')
                                    length_C = g.nodes[C]['length']
                                else:
                                    length_C = sexSize(path_data, gene,
                                                       C)[0][3]
                                if length_C <= 5:  #si le candidat est de taille inf a 5 je l'ajoute
                                    print('C IS A NODE, INF TO 5',
                                          g.has_node(C))
                                    print('extension inf a 5', eve, pair, X, C,
                                          direction)
                                    #flag=False
                                    isgood = codage(eve2pair, C, d_X)
                                    if isgood:
                                        path = '/'.join(path)
                                        #on update g.nodes()
                                        if direction == 'Cter':  #and (X==A or X==B) interet du ET ?
                                            update = X + '$.' + C
                                            path = path.replace(
                                                f'{X}/{C}', update)
                                        elif direction == 'Nter':
                                            update = C + '.' + X + '$'
                                            path = path.replace(
                                                f'{C}/{X}', update)
                                        eve2pair[eve][-1][path_bool] = path
                                        path = path.split('/')
                                        mem_nodes[X] = update
                                        mapping = {X: update}
                                        print('MAPPING', mapping)
                                        g = nx.relabel_nodes(
                                            g, mapping
                                        )  #j'update le noeud seed X du graphe en le noeud seed-candidat
                                        g.nodes[update]['length'] += length_C
                                        marge -= length_C
                                        print('after', marge)
                                        try:  # l'unite du type { X , C , ...} devient { X/C , C , ...} après extension, il reste à enlever le doublon C
                                            g.remove_node(C)
                                            mem_nodes[C] = update
                                        except nx.exception.NetworkXError:  # cas où C a déjà été enlevé
                                            pass
                                    else:
                                        break
                                else:  #si le noeud est de taille >= 5 je dois regarder calculer les hhr de X et celui en face de X et du candidat avec celui en face de X
                                    #G_init=sim_graph_workframe(gene)
                                    taille_A = g.nodes[A]['length']
                                    colA_bis, colB_bis = get_alnCols(
                                        g, A, B, eveDup)
                                    print('ZZZZZZZZZZZBEFORE', A, B, colA_bis,
                                          colB_bis)
                                    if X == A:  #dans la paire (A,B) la seed X est A
                                        if B != C:
                                            if g.has_edge(C, B):
                                                colB = get_alnCols(
                                                    g, C, B, eveDup)[1]
                                            else:
                                                break
                                            if (float(colB_bis[0]) <= float(
                                                    colB[0]) <= float(
                                                        colB_bis[1])
                                                    or float(colB_bis[0]) <=
                                                    float(colB[1]) <=
                                                    float(colB_bis[1])) or (
                                                        float(colB[0]) <=
                                                        float(colB_bis[0]) <=
                                                        float(colB[1])
                                                        or float(colB[0]) <=
                                                        float(colB_bis[1]) <=
                                                        float(colB[1])):
                                                break  #on verifie que lalignement entre A et B et l'alignement correspondant impliquant l'extension ne s'overlap pas
                                            else:
                                                isgood = codage(
                                                    eve2pair, C, d_X)
                                                if isgood:
                                                    path = '/'.join(path)
                                                    if direction == 'Cter':
                                                        if float(
                                                                colB_bis[1]
                                                        ) < float(
                                                                colB[0]
                                                        ):  #on verifie que les bornes de lalignement sont coherentes
                                                            update = X + '$.' + C
                                                            path = path.replace(
                                                                f'{X}/{C}',
                                                                update)
                                                        else:
                                                            break

                                                    elif direction == 'Nter':
                                                        if float(
                                                                colB[1]
                                                        ) < float(colB_bis[0]):
                                                            update = C + '.' + X + '$'
                                                            path = path.replace(
                                                                f'{C}/{X}',
                                                                update)
                                                        else:
                                                            break
                                                else:
                                                    break

                                                for edge in g.edges(X):
                                                    if direction == 'Cter':
                                                        temp = g.edges[edge][
                                                            'info'][-1].split(
                                                                '-')
                                                        print(
                                                            'ZZZZZZZZZZZ',
                                                            temp, edge, X)
                                                        temp[1] = str(
                                                            int(temp[1]) +
                                                            length_C)
                                                        print(
                                                            'ZZZZZZZZZZZOA',
                                                            temp, edge, X)
                                                        g.edges[edge]['info'][
                                                            -1] = '-'.join(
                                                                temp)
                                                    if direction == 'Nter':
                                                        temp = g.edges[edge][
                                                            'info'][-1].split(
                                                                '-')
                                                        print(
                                                            'ZZZZZZZZZZZ',
                                                            temp)
                                                        temp[0] = str(
                                                            int(temp[0]) -
                                                            length_C)
                                                        if int(
                                                                temp[0]
                                                        ) - length_C < 0:
                                                            temp[0] = str(1)
                                                        print(
                                                            'ZZZZZZZZZZZOA',
                                                            temp, edge, X)
                                                        g.edges[edge]['info'][
                                                            -1] = '-'.join(
                                                                temp)
                                                eve2pair[eve][-1][
                                                    path_bool] = path
                                                path = path.split('/')
                                                mem_nodes[X] = update
                                                mapping = {
                                                    X: update,
                                                    C: update
                                                }
                                                print('MAPPING', mapping)
                                                g = nx.relabel_nodes(
                                                    g, mapping)
                                                g.nodes[update][
                                                    'length'] += length_C
                                                marge -= length_C
                                                print('after', marge)
                                                try:
                                                    #mapping={C:update}
                                                    mem_nodes[C] = update
                                                    #g=nx.relabel_nodes(g, mapping)
                                                    #g.remove_node(C)
                                                except KeyError:
                                                    pass

                                        else:
                                            break
                                    else:  #dans la paire (A,B) la seed X est B
                                        if A != C:
                                            if g.has_edge(A, C):
                                                colA = get_alnCols(
                                                    g, A, C, eveDup)[0]
                                            else:
                                                print('AAAAA')
                                                print('mem_nodes', mem_nodes)
                                                break
                                            if (float(colA_bis[0]) <= float(
                                                    colA[0]) <= float(
                                                        colA_bis[1])
                                                    or float(colA_bis[0]) <=
                                                    float(colA[1]) <=
                                                    float(colA_bis[1])) or (
                                                        float(colA[0]) <=
                                                        float(colA_bis[0]) <=
                                                        float(colA[1])
                                                        or float(colA[0]) <=
                                                        float(colA_bis[1]) <=
                                                        float(colA[1])):
                                                break
                                            else:
                                                isgood = codage(
                                                    eve2pair, C, d_X)
                                                path = '/'.join(path)
                                                if isgood:
                                                    if direction == 'Cter':
                                                        if float(colA_bis[1]
                                                                 ) < float(
                                                                     colA[0]):
                                                            update = X + '$.' + C  #ITOU WARNING
                                                            path = path.replace(
                                                                f'{X}/{C}',
                                                                update)
                                                        else:
                                                            break

                                                    if direction == 'Nter':
                                                        if float(
                                                                colA[1]
                                                        ) < float(colA_bis[0]):
                                                            update = C + '.' + X + '$'  #ATTENTION ICI AVANT A CT B
                                                            path = path.replace(
                                                                f'{C}/{X}',
                                                                update)
                                                        else:
                                                            break

                                                else:
                                                    break

                                                for edge in g.edges(X):
                                                    if direction == 'Cter':
                                                        temp = g.edges[edge][
                                                            'info'][-2].split(
                                                                '-')
                                                        print(
                                                            'ZZZZZZZZZZZ',
                                                            temp, edge, X)
                                                        temp[1] = str(
                                                            int(temp[1]) +
                                                            length_C)
                                                        print(
                                                            'ZZZZZZZZZZZO',
                                                            temp, edge, X)
                                                        g.edges[edge]['info'][
                                                            -1] = '-'.join(
                                                                temp)
                                                    if direction == 'Nter':
                                                        temp = g.edges[edge][
                                                            'info'][-2].split(
                                                                '-')
                                                        print(
                                                            'ZZZZZZZZZZZ',
                                                            temp, edge, X)
                                                        temp[0] = str(
                                                            int(temp[0]) -
                                                            length_C)
                                                        if int(
                                                                temp[0]
                                                        ) - length_C < 0:
                                                            temp[0] = str(1)
                                                        print(
                                                            'ZZZZZZZZZZZO',
                                                            temp, edge, X)
                                                        g.edges[edge]['info'][
                                                            -1] = '-'.join(
                                                                temp)
                                                eve2pair[eve][-1][
                                                    path_bool] = path
                                                path = path.split('/')
                                                mem_nodes[X] = update
                                                mapping = {
                                                    X: update,
                                                    C: update
                                                }
                                                print('MAPPING', mapping)
                                                g = nx.relabel_nodes(
                                                    g, mapping)
                                                g.nodes[update][
                                                    'length'] += length_C
                                                marge -= length_C
                                                print('after', marge)
                                                try:
                                                    #mapping={C:update}
                                                    mem_nodes[C] = update
                                                    #g=nx.relabel_nodes(g, mapping)
                                                    #g.remove_node(C)
                                                except KeyError:
                                                    pass

                                        else:
                                            break

                    else:
                        break

                    print('X end while before', X)
                    X = update
                    print('X end while after', X)
                    flag = True
                    print('end while, g.hasnodeA/B', g.has_node(A),
                          g.has_node(A))
                    if not g.has_node(A):
                        A = mem_nodes[A]
                    if not g.has_node(B):
                        B = mem_nodes[B]
                    print('endwhile A/B mem_nodes update', A, B)

                lst = [[key, mem_nodes[key]] for key in mem_nodes]
                lst_graph = nx.from_edgelist(lst, create_using=nx.DiGraph)
                for node in lst_graph.nodes:
                    if lst_graph.out_degree(node) == 0:
                        for node_bis in lst_graph.nodes:
                            if node_bis in node and node_bis != node:
                                mem_nodes[node_bis] = node

    return g


def writecsv_ASRU(path_data, gene):
    with open(f'{path_data}/data/{gene}/{gene}_ASRUs_table.csv',
              'w') as f, open(
                  f'{path_data}/data/{gene}/{gene}_instances_table.csv',
                  'w') as g:
        csv_writer1 = csv.writer(f, delimiter=',')
        csv_writer2 = csv.writer(g, delimiter=',')

        csv_writer1.writerow([
            'gene', 'ASRU', 'Nbinstances', 'max', 'min', 'moy', 'median',
            'std', 'eventsRank'
        ])
        csv_writer2.writerow(['instance', 'size', 'NbSex', 'ASRU', 'gene'])

        G_sim_gene = extensionbis(path_data, gene)
        #G_original_gene=sim_graph_initial_bis(path_data, gene)
        ConnectedComp_gene = nx.connected_components(G_sim_gene)
        for ASRU in ConnectedComp_gene:
            gene_csv_line1 = []
            gene_csv_line1.append(gene)  #nom du gene
            gene_csv_line1.append(ASRU)
            gene_csv_line1.append(len(ASRU))
            instancedatatmp = []
            eve = []
            for pair in list(itertools.combinations(list(ASRU), 2)):
                if (pair[0], pair[1]) in G_sim_gene.edges():
                    eve += G_sim_gene.edges[pair]['info'][0]
            eve = list(np.unique(np.array(eve)))
            for inst in ASRU:
                gene_csv_line2 = []
                gene_csv_line2.append(inst)
                if '.' not in inst:
                    instancedatatmp.append(
                        sexSize(path_data, gene, inst)[0][3])
                    gene_csv_line2.append(sexSize(path_data, gene, inst)[0][3])
                    gene_csv_line2.append(1)
                else:
                    sizesum = 0
                    for seed in inst.split('.'):
                        if '$' not in seed:
                            sizesum += sexSize(path_data, gene, seed)[0][3]
                        else:
                            sizesum += sexSize(path_data, gene,
                                               seed.rstrip('$'))[0][3]
                    gene_csv_line2.append(sizesum)
                    gene_csv_line2.append(len(inst.split('.')))
                    instancedatatmp.append(sizesum)
                gene_csv_line2.append(ASRU)
                gene_csv_line2.append(gene)
                csv_writer2.writerow(gene_csv_line2)
            gene_csv_line1.append(max(instancedatatmp))
            gene_csv_line1.append(min(instancedatatmp))
            gene_csv_line1.append(np.mean(instancedatatmp))
            gene_csv_line1.append(np.median(instancedatatmp))
            gene_csv_line1.append(np.std(instancedatatmp))
            gene_csv_line1.append(eve)
            csv_writer1.writerow(gene_csv_line1)


def run():
    args = parse_args(sys.argv[1:])
    gene = args.geneName
    path_data = args.dataPATH
    writecsv_ASRU(path_data, gene)


if __name__ == '__main__':
    run()