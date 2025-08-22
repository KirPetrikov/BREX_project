import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import os
import csv


##################################
# Counting how many hits DF&PADLOC#
##################################

def parse_padloc_and_defensefinder_output(accession, output_dir, padloc_dir, df_dir,
                                          conflict_winner, df_err, pl_err,
                                          ref_file=Path(os.path.dirname(sys.argv[0]),
                                                        "immune_system_list_reference.csv"),
                                          warning_file="warnings.txt"):
    warning_file = os.path.join(output_dir, warning_file)
    ref_dict = load_reference_file(ref_file)

    padloc_db_dict = parse_padloc(accession, padloc_dir, ref_dict, pl_err)

    df_db_dict = parse_df(accession, df_dir, ref_dict, df_err)  # ПРЕДПОЛАГАЕТСЯ, ЧТО ACCESSION - ОДИН ЭЛЕМЕНТ

    combine_df_padloc_output(accession, padloc_db_dict, df_db_dict, output_dir, warning_file, conflict_winner)


def load_reference_file(ref_file):
    # Load file with all the name differences between PADLOC and Defence Finder
    with open(ref_file) as ref_list_file:
        ref_list = csv.reader(ref_list_file)
        start = True
        shared_systems = []  # list with systems both softwares look for
        ref_dict = {}  # dictionary with all reference values to equalise the systems
        for line in ref_list:
            if not start:
                if "" not in line:
                    shared_systems.append(line[0])
                if line[1] != "":
                    if "/" not in line[1]:  # multiple names
                        ref_dict[line[1]] = line[0]
                    else:
                        names = line[1].split("/")
                        for name in names:
                            ref_dict[name] = line[0]
                if line[2] != "":
                    if "/" not in line[2]:  # multiple names
                        ref_dict[line[2]] = line[0]
                    else:
                        names = line[2].split("/")
                        for name in names:
                            ref_dict[name] = line[0]
            start = False
    return ref_dict


def parse_padloc(accession, padloc_dir, ref_dict, pl_err):
    ###parse PADLOC

    padloc_db_dict = {}
    for padloc_file in os.listdir(padloc_dir):
        if padloc_file.endswith(".csv"):
            padloc_db_dict[accession] = []
            padloc_file = open(os.path.join(padloc_dir, padloc_file), mode="r", encoding='utf-8-sig')
            padloc_list = csv.reader(padloc_file, quotechar='"', delimiter=",", quoting=csv.QUOTE_ALL,
                                     skipinitialspace=True)

            start = True
            system_numbers = []
            systems_dict = {}
            for line_list in padloc_list:
                if start:
                    for i in range(len(line_list)):
                        if line_list[i] == "system":
                            system_index = i
                        if line_list[i] == "system.number":
                            system_number_index = i
                        if line_list[i] == "target.name":
                            gene_index = i
                        if line_list[i] == "protein.name":
                            protein_name_index = i
                    start = False
                else:
                    # get relevant values
                    system = line_list[system_index].split("_")[0]
                    subtype = line_list[system_index]
                    gene = line_list[gene_index]
                    protein_name = line_list[protein_name_index]

                    # check if this instance of the system was already added
                    system_number = line_list[system_number_index]
                    if subtype != "DRT_other":  # these systems also exist as retron XIII and are identified solely as such by the online padloc version
                        if system_number not in systems_dict:
                            if system in ref_dict:
                                system = ref_dict[system]
                            else:
                                pl_err.add(system)
                                # print ("padloc ref error:", system)
                            if system != "skip":  # Like DMS with look like these are unfinished potential partial systems (or just plain duplicates), not necessary to investigate on full genomes like these
                                systems_dict[system_number] = [system, subtype, [gene], [protein_name]]
                        else:
                            systems_dict[system_number][2] = systems_dict[system_number][2] + [gene]
                            systems_dict[system_number][3] = systems_dict[system_number][3] + [protein_name]
            # remove systems identified more than once that are identified as the same master family (happens for e.g. nhi and AbiO which are counted as the same system by DefenseFinder)
            temp = []
            res = dict()
            for key, val in systems_dict.items():
                if val not in temp:
                    temp.append(val)
                    res[key] = val
            systems_dict = res.copy()
            # make final list of systems
            for key in systems_dict:
                system = systems_dict[key][0]
                subtype = systems_dict[key][1]
                gene_list = ";".join(systems_dict[key][2])
                protein_name_list = ";".join(systems_dict[key][3])
                padloc_db_dict[accession] = padloc_db_dict[accession] + [
                    [system, subtype, gene_list, protein_name_list]]
            padloc_file.close()
    return padloc_db_dict


def parse_df(accession, df_dir, ref_dict, df_err):
    # parse defence finder

    df_db_dict = {}
    for df_file in os.listdir(df_dir):
        if df_file[-len('_systems.tsv'):] == '_systems.tsv':
            df_db_dict[accession] = []
            df = open(os.path.join(df_dir, df_file), mode="r", encoding='utf-8-sig')
            df_list = csv.reader(df, quotechar='"', delimiter="\t", quoting=csv.QUOTE_ALL, skipinitialspace=True)

            start = True
            for line_list in df_list:
                if start:
                    for i in range(len(line_list)):
                        if line_list[i] == "type":
                            system_index = i
                        elif line_list[i] == "subtype":
                            subtype_index = i
                        elif line_list[i] == "protein_in_syst":
                            genes_index = i
                        elif line_list[i] == "name_of_profiles_in_sys":
                            gene_names_index = i

                    start = False
                else:
                    # get relevant values
                    system = line_list[system_index]
                    if system in ref_dict:
                        if system != "skip":
                            system = ref_dict[system]
                    else:
                        df_err.add(system)
                        # print ("defence finder ref error:", system)
                    subtype = line_list[subtype_index]
                    gene_list = ";".join(line_list[genes_index].split(","))
                    gene_names = line_list[gene_names_index]

                    df_db_dict[accession] = df_db_dict[accession] + [[system, subtype, gene_list, gene_names]]
    return df_db_dict


def combine_df_padloc_output(accession, padloc_db_dict, df_db_dict, output_dir,
                             warning_file, conflict_winner):
    combined_dict = {}

    PL_dict = {}
    if accession in padloc_db_dict.keys():
        for system in padloc_db_dict[accession]:
            PL_dict[";".join(system)] = system

    DF_dict = {}
    if accession in df_db_dict.keys():
        for system in df_db_dict[accession]:
            DF_dict[";".join(system)] = system

    PL_dict = remove_internal_conflict(PL_dict, DF_dict, "PL")
    DF_dict = remove_internal_conflict(DF_dict, PL_dict, "DF")

    if accession in padloc_db_dict.keys():
        for system in PL_dict.keys():
            system_db = PL_dict[system]
            combined_dict = test_system_for_overlap(system_db, combined_dict, "PL", DF_dict, conflict_winner, accession,
                                                    warning_file)

    if accession in df_db_dict.keys():
        for system in DF_dict.keys():
            system_db = DF_dict[system]
            combined_dict = test_system_for_overlap(system_db, combined_dict, "DF", PL_dict, conflict_winner, accession,
                                                    warning_file)

    with open(Path(output_dir, accession + ".csv"), "w") as output_file:
        output_file.write(
            "system family DefenseFinder,system family PADLOC,subtype DefenseFinder,subtype PADLOC,genes, genes DefenseFinder, genes PADLOC\n")

        # combine data
        for genes in combined_dict.keys():
            info = combined_dict[genes]

            # make summary file for accession
            system_PL = "N.A."
            system_DF = "N.A."
            subtype_PL = "N.A."
            subtype_DF = "N.A."
            genes_DF = "N.A."
            genes_PL = "N.A."

            if "DF" in info.keys():
                system_DF = info["DF"][0]
                subtype_DF = info["DF"][1]
                genes_DF = info["DF"][3].replace(",", ";")

            if "PL" in info.keys():
                system_PL = info["PL"][0]
                subtype_PL = info["PL"][1]
                genes_PL = "N.A."
                genes_PL = info["PL"][3].replace(",", ";")

            output_file.write(
                "%s,%s,%s,%s,%s,%s,%s\n" % (system_DF, system_PL, subtype_DF, subtype_PL, genes, genes_DF, genes_PL))

    return


def remove_internal_conflict(db_from_self, db_from_other, testing):
    new_db_from_self = dict()

    # test internal conflict (systems with differen families assigned but the exact same genes within the DF or PL database)
    removed_list = []
    flag = True
    for system1 in db_from_self.keys():
        remove = False
        internal_conflict = False
        for system2 in db_from_self.keys():
            genes_1 = db_from_self[system1][2]
            genes_2 = db_from_self[system2][2]
            genes_set_1 = set(genes_1.split(";"))
            genes_set_2 = set(genes_2.split(";"))
            if genes_set_1 == genes_set_2:
                if system1 != system2:
                    internal_conflict = True
                    db = db_from_self[system1]
                    accession = "not_main"
                    warning_file = "not_main"
                    found_by_both_b, part_of_a_whole, part_of_a_whole_longest, overlap_not_part_of_a_whole, different_family_same_system, overlap_not_part_of_a_whole, different_family, no_conflict = find_overlap(
                        db, genes_set_1, db_from_other, accession, testing, warning_file)
                    if genes_set_1 in removed_list:  # previously merged
                        remove = True
                    elif db[0] == db_from_self[system2][0]:  # same system family, differen subsystem: merge
                        new_db_from_self[system1] = db[0], db[1] + "/" + db_from_self[system2][1], db[2], db[3] + "/" + \
                                                           db_from_self[system2][3]
                        removed_list.append(genes_set_1)
                    elif found_by_both_b:
                        new_db_from_self[system1] = db_from_self[system1]  # keep unchanged
                    elif different_family_same_system:
                        remove = True  # do not test
                    else:  # unique find, merge
                        new_db_from_self[system1] = db[0] + "/" + db_from_self[system2][0], db[1] + "/" + \
                                                    db_from_self[system2][1], db[2], db[3] + "/" + \
                                                    db_from_self[system2][3]
                        removed_list.append(genes_set_1)
        if internal_conflict == False:
            new_db_from_self[system1] = db_from_self[system1]
    return new_db_from_self


def find_overlap(db, genes_set, db_from_other, accession, testing, warning_file):
    no_conflict = False
    found_by_both_b = False
    part_of_a_whole_longest = False
    part_of_a_whole = False
    overlap_not_part_of_a_whole = False
    different_family_same_system = False
    different_family = False
    internal_conflict = False
    overlap = False

    # test for conflict between PL and DF
    conflicting_systems = []
    for system in db_from_other.keys():
        genes_from_other = db_from_other[system][2]
        genes_from_other_set = set(genes_from_other.split(";"))
        combined_genes_set = genes_set.union(genes_from_other_set)
        if len(combined_genes_set) < (len(genes_set) + len(genes_from_other_set)):
            overlap = True
            conflicting_systems.append(db_from_other[system])  # define which system(s) conflict

    if overlap:  # Does the system have overlap with a system from the other?
        for conflicting_system in conflicting_systems:
            genes_from_other_set = set(conflicting_system[2].split(";"))
            combined_genes_set = genes_set.union(genes_from_other_set)
            if db[0] == conflicting_system[0]:  # Do these systems share the same system family?
                if genes_set == genes_from_other_set:  # Do they have exactly the same genes?
                    found_by_both_b = True
                else:
                    if (combined_genes_set == genes_set) or (
                            combined_genes_set == genes_from_other_set):  # is one contained in the other?
                        part_of_a_whole = True
                        if combined_genes_set == genes_set:  # keep longest genes
                            part_of_a_whole_longest = True
                    else:
                        overlap_not_part_of_a_whole = True

            else:
                if genes_set == genes_from_other_set:  # do they have exactly the same genes?
                    different_family_same_system = True
                    with open(warning_file, "a+") as w_f:
                        w_f.write(
                            "%s,%s,%s,%s\n" % (
                            accession, db[0], conflicting_system[0], db[2]))
                else:
                    overlap_not_part_of_a_whole = True
                    different_family = True
    else:
        no_conflict = True

    return found_by_both_b, part_of_a_whole, part_of_a_whole_longest, overlap_not_part_of_a_whole, different_family_same_system, overlap_not_part_of_a_whole, different_family, no_conflict


def test_system_for_overlap(db, combined_dict, testing, db_from_other, conflict_winner, accession, warning_file):
    genes = db[2]
    genes_set = set(genes.split(";"))
    conflicting_systems = []

    found_by_both_b, part_of_a_whole, part_of_a_whole_longest, overlap_not_part_of_a_whole, different_family_same_system, overlap_not_part_of_a_whole, different_family, no_conflict = find_overlap(
        db, genes_set, db_from_other, accession, testing, warning_file)

    tester_wins = False
    conflict_winner_wins = False
    if found_by_both_b:
        tester_wins = True
    elif no_conflict:
        tester_wins = True
    elif part_of_a_whole:
        if part_of_a_whole_longest:
            tester_wins = True
    elif different_family_same_system or overlap_not_part_of_a_whole:
        conflict_winner_wins = True

    if conflict_winner_wins:
        if testing == "PL" and conflict_winner == "PL":
            if genes in combined_dict.keys():
                combined_dict[genes]["PL"] = db
            else:
                combined_dict[genes] = {"PL": db}
        elif testing == "DF" and conflict_winner == "DF":
            if genes in combined_dict.keys():
                combined_dict[genes]["DF"] = db
            else:
                combined_dict[genes] = {"DF": db}
    elif tester_wins:
        if testing == "PL":
            if genes in combined_dict.keys():
                combined_dict[genes]["PL"] = db
            else:
                combined_dict[genes] = {"PL": db}
        else:
            if genes in combined_dict.keys():
                combined_dict[genes]["DF"] = db
            else:
                combined_dict[genes] = {"DF": db}
    return combined_dict


def find_immune_systems(sample_id, df_err,
                        pl_err,
                        df_dir=None,
                        pl_dir=None,
                        output_dir=None,
                        conflict_winner='DF',
                        cpu=None):
    #############
    # parse input#
    #############

    accession = sample_id  # ДЛЯ ОДНОГО ОБРАЗЦА

    if not output_dir:
        output_dir = Path(Path.cwd(), "PDF_comb_results")
        output_dir.mkdir(parents=True, exist_ok=True)
    if not output_dir.is_dir():
        raise FileNotFoundError(f'No such directory: {output_dir}')

    if (conflict_winner != "DF") and (conflict_winner != "PL"):
        raise Exception("Conflict winner not recognised. Please provide DF or PL.")

    ##############
    # parse output#
    ##############
    parse_padloc_and_defensefinder_output(accession, output_dir, pl_dir, df_dir,
                                          conflict_winner, df_err, pl_err)


def step1_one_sample(idx, df_pl_viz):
    df_err = set()  # ожидается, что оно останется пустым
    pl_err = set()  # ожидается, что оно останется пустым

    df_dir = df_pl_viz.loc[idx, 'df_dir']
    if not Path(df_dir).exists():
        raise Exception(f'Path {df_dir} not found')

    pl_dir = df_pl_viz.loc[idx, 'pl_dir']
    if not Path(pl_dir).exists():
        raise Exception(f'Path {pl_dir} not found')

    find_immune_systems(str(df_pl_viz.loc[idx, 'sample_id']), df_err, pl_err, df_dir=df_dir, pl_dir=pl_dir)

    if df_err:
        print('Defensefinder reference ERROR, please update immune_system_list_reference.csv:',
              str(df_pl_viz.loc[idx, 'sample_id']),
              df_err)
    if pl_err:
        print('Padloc reference ERROR, please update immune_system_list_reference.csv:',
              str(df_pl_viz.loc[idx, 'sample_id']),
              pl_err)
    return


def step1(df_pl_viz):
    num_threads = 1  # Оптимальное количество потоков

    exceptions_socket = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        future_to_i = {executor.submit(step1_one_sample, idx, df_pl_viz): idx for idx in df_pl_viz.index}
        for future in as_completed(future_to_i):
            i = future_to_i[future]
            try:
                future.result()
            except Exception as exc:
                print(f'Образец {i} вызвал исключение: {exc}')
                exceptions_socket.append([i, exc])
                # with open('results/exceptions_step1.txt', 'w') as file_ex:
                #     file_ex.write('\n'.join([''.join(l2) for l2 in exceptions_socket]))
    print('FINISHED')
    if exceptions_socket:
        print(f'Образцы, вызвавшие исключения: {exceptions_socket}')
    return
