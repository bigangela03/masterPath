package masterPATH;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Wrapper class contains wrappers implementing the whole pipeline
 *
 * @author Natalia Rubanova
 */
public class Wrapper {

    DBManager dbutils = new DBManager();
    NetworkManager nutils = new NetworkManager();
    PathwayManager putils = new PathwayManager();
    CentralityManager outils = new CentralityManager();
    NetworkTopology utils = new NetworkTopology();
    Wrapper wrap = new Wrapper();
    RandomManager rand = new RandomManager();

    public Wrapper() throws IOException {

    }

    public Network load_network(String file_nodes, String file_interaction) throws IOException {
        Network nw = new Network();
        nw.loadNetworkfromfile(file_nodes, file_interaction);
        nutils.build_map_of_neighbors(nw);
        return nw;
    }


    public void find_shortest_paths_and_calculate_centrality(
            Network nw,
            String file_hitlist,
            String file_fimplementers,
            String file_output,
            String prefix,
            int max_length_for_shortest_path,
            int min_len_for_path,
            int max_len_for_path,
            String folder_for_random_paths,
            String prefix_for_random_paths,
            int number_of_permutations) throws IOException {

        nutils.load_hitlist_and_finalimpl(file_hitlist, file_fimplementers, dbutils.hugo);

        nutils.find_pathway_for_list_BF_algorithm(nw, dbutils.hugo,
                max_length_for_shortest_path, file_output, prefix);

        putils.find_the_shortest_paths(file_output, "_shortest_paths",
                nw, nutils.hg, nutils.fpl, 0, 0, 0, 0);
        putils.add_connectivity_to_pathways(nw, file_output + "_shortest_paths",
                file_output + "_shortest_paths_with_connectivity_tmp");

        //
        outils.calculate_centrality_scores_for_paths(nw,
                file_output + "_shortest_paths",
                file_output + "_paths_centrality_tmp",
                min_len_for_path, max_len_for_path);
        outils.add_hgnc_symbols_to_paths(file_output + "_paths_centrality_tmp",
                file_output + "_paths_centrality_w_hugo_tmp",
                nw, dbutils.hugo_by_id);
        outils.calculate_paths_degree(file_output + "_paths_centrality_w_hugo_tmp",
                file_output + "_shortest_paths_with_connectivity_tmp",
                file_output + "_paths_centrality");

        outils.calculate_centrality_scores_for_nodes(nw, dbutils.hugo_by_id,
                nutils, file_output + "_shortest_paths",
                file_output + "_nodes_centrality");

        File file = new File(file_output + "_shortest_paths_with_connectivity_tmp");
        file.delete();
        file = new File(file_output + "_paths_centrality_tmp");
        file.delete();
        file = new File(file_output + "_paths_centrality_w_hugo_tmp");
        file.delete();

        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list = new ArrayList(nutils.hg.keySet());

        rand.permute_phenotype_label(nw, nutils, putils, outils, dbutils,
                folder_for_random_paths, prefix_for_random_paths, number_of_permutations,
                hit_list, dbutils.hugo_by_id, file_fimplementers);

        rand.calculate_p_values_phenotype_label_permutation(nw, nutils, putils,
                outils, dbutils,
                folder_for_random_paths + "_random_paths_shortest",
                file_output + "_paths_centrality",
                folder_for_random_paths + "random_hit_lists",
                initial_hitlist_size, max_len_for_path + 2, dbutils.hugo_by_id);
    }

    public void find_shortest_paths_and_calculate_centrality_ppi(
            Network nw,
            String file_hitlist,
            String file_fimplementers,
            String file_output,
            String prefix,
            int max_length_for_shortest_path,
            int min_len_for_path,
            int max_len_for_path,
            String folder_for_random_paths,
            String prefix_for_random_paths,
            int number_of_permutations) throws IOException {

        nutils.load_hitlist_and_finalimpl(file_hitlist, file_fimplementers, dbutils.hugo);

        nutils.find_pathway_for_list_BF_algorithm_ppi(nw, dbutils.hugo,
                max_length_for_shortest_path, file_output, prefix);

        putils.find_the_shortest_paths(file_output, "_shortest_paths",
                nw, nutils.hg, nutils.fpl, 0, 0, 0, 0);
        putils.add_connectivity_to_pathways(nw, file_output + "_shortest_paths",
                file_output + "_shortest_paths_with_connectivity_tmp");

        //
        outils.calculate_centrality_scores_for_paths_ppi(nw,
                file_output + "_shortest_paths",
                file_output + "_paths_centrality_tmp",
                min_len_for_path, max_len_for_path);
        outils.add_hgnc_symbols_to_paths(file_output + "_paths_centrality_tmp",
                file_output + "_paths_centrality_w_hugo_tmp",
                nw, dbutils.hugo_by_id);
        outils.calculate_paths_degree(file_output + "_paths_centrality_w_hugo_tmp",
                file_output + "_shortest_paths_with_connectivity_tmp",
                file_output + "_paths_centrality");

        outils.calculate_centrality_scores_for_nodes(nw, dbutils.hugo_by_id,
                nutils, file_output + "_shortest_paths",
                file_output + "_nodes_centrality");

        File file = new File(file_output + "_shortest_paths_with_connectivity_tmp");
        file.delete();
        file = new File(file_output + "_paths_centrality_tmp");
        file.delete();
        file = new File(file_output + "_paths_centrality_w_hugo_tmp");
        file.delete();

        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list = new ArrayList(nutils.hg.keySet());

        rand.permute_phenotype_label(nw, nutils, putils, outils, dbutils,
                folder_for_random_paths, prefix_for_random_paths, number_of_permutations,
                hit_list, dbutils.hugo_by_id, file_fimplementers);

        rand.calculate_p_values_phenotype_label_permutation_ppi(nw, nutils, putils,
                outils, dbutils,
                folder_for_random_paths + "_random_paths_shortest",
                file_output + "_paths_centrality",
                folder_for_random_paths + "random_hit_lists",
                initial_hitlist_size, max_len_for_path + 2, dbutils.hugo_by_id);

        rand.calculate_p_values_phenotype_label_permutation_nodes(nw, nutils, putils,
                outils, dbutils,
                folder_for_random_paths + "_random_paths_shortest",
                file_output + "_nodes_centrality",
                folder_for_random_paths + "random_hit_lists",
                initial_hitlist_size,
                dbutils.hugo_by_id);

    }


    /**
     * Part 1 of pipeline for DNA repair screen in ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void ogg1_part1_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlOG, FoldersPaths.fpOG_ppi, dbutils.hugo);
        nutils.find_pathway_for_list_BF_algorithm_ppi(all, dbutils.hugo, 5, FoldersPaths.foundDFOG_ppi, "OG");
        putils.find_the_shortest_paths(FoldersPaths.foundDFOG_ppi, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     * Part 2 of pipeline for DNA repair screen in ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void ogg1_part2_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.add_connectivity_to_pathways(all, FoldersPaths.foundDFOG_ppi + "_top_ranked", FoldersPaths.foundDFOG_ppi + "_top_ranked_with_connectivity");
        //
        outils.calculate_centrality_scores_for_paths_ppi(all, FoldersPaths.foundDFOG_ppi + "_top_ranked", FoldersPaths.foundDFOG_ppi + "_top_ranked__unique_overrepr", 2, 3);
        outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFOG_ppi + "_top_ranked__unique_overrepr", FoldersPaths.foundDFOG_ppi + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.calculate_paths_degree(FoldersPaths.foundDFOG_ppi + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFOG_ppi + "_top_ranked_with_connectivity", FoldersPaths.foundDFOG_ppi + "_top_ranked__unique_overrepr_w_names_pthconn");

    }

    /**
     * Part 'p-values' of pipeline for DNA repair screen in ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void ogg1_random_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlOG, FoldersPaths.fpOG_ppi, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list = new ArrayList(nutils.hg.keySet());
        rand.permute_phenotype_label(all, nutils, putils, outils, dbutils, FoldersPaths.OGrandom_path_ppi, "OGR", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpOG_ppi);
        rand.calculate_p_values_phenotype_label_permutation_ppi(all, nutils, putils, outils, dbutils, FoldersPaths.OGrandom_path_ppi + "foundDF.txt_top_ranked", FoldersPaths.foundDFOG_ppi + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.OGrandom_path_ppi + "RandomHitGenes_lists.txt", initial_hitlist_size, 5, dbutils.hugo_by_id);

    }

    /**
     * Part 4 of pipeline for DNA repair screen in ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void ogg1_part4_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlOG, FoldersPaths.fpOG_ppi, dbutils.hugo);
        putils.create_cyto_for_pathways(FoldersPaths.foundDFOG_ppi + "_top_ranked", FoldersPaths.outputFOG + "ppi\\networks\\paths.txt", FoldersPaths.outputFOG + "ppi\\networks\\", all, nutils.hg, nutils.fpl);
    }

    /**
     * Part 1 of pipeline for miRNA muscle differentiation screen in mixed
     * network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void mirna63Sys_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        nutils.find_pathway_for_list_BF_algorithm(all, dbutils.hugo, 6, FoldersPaths.foundDF63, "MR");
        putils.find_the_shortest_paths(FoldersPaths.foundDF63, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     * Part 2 of pipeline for miRNA muscle differentiation screen in mixed
     * network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void mirna63Sys_part2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.add_connectivity_to_pathways(all, FoldersPaths.foundDF63 + "_top_ranked", FoldersPaths.foundDF63 + "_top_ranked_with_connectivity");
        //
        outils.calculate_centrality_scores_for_paths(all, FoldersPaths.foundDF63 + "_top_ranked", FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr", 1, 4);
        outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr", FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.calculate_paths_degree(FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDF63 + "_top_ranked_with_connectivity", FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr_w_names_pthconn");

        outils.calculate_centrality_scores_for_nodes(all, dbutils.hugo_by_id, nutils, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked_hubs");

    }

    /**
     * Part 'p-values' of pipeline for miRNA muscle differentiation screen in
     * mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void mirna63Sys_random(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list = new ArrayList(nutils.hg.values());

        for (String s : nutils.hg.keySet()) {
            //System.out.println(s + " " + nutils.hg.get(s));
        }
        //rand.check_distributions(all, dbutils.hugo_by_id);
        //rand.permute_phenotype_label(all, nutils, putils, outils, dbutils, FoldersPaths.M63random_path, "63R", 2000, hit_list, dbutils.hugo_by_id, FoldersPaths.fp63);
        //rand.calculate_p_values_phenotype_label_permutation(all, nutils, putils, outils, dbutils, FoldersPaths.M63random_path + "foundDF.txt_top_ranked", FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.M63random_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);
        String fn = "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\names_conversion.txt";
        //outils.add_experimnetal_scores_and_aggregated_pvalues_to_paths(FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr_w_names_pthconn"+ "_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\expression.txt", fn, FoldersPaths.foundDF63 + "_top_ranked_aggr_pvalues");
        //outils.add_experimnetal_scores_and_aggregated_pvalues_to_paths(FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names_pthconn"+ "_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\expression.txt", fn, FoldersPaths.foundDFC2 + "_top_ranked_aggr_pvalues");
        outils.add_experimental_scores_and_pvalues_to_nodes(FoldersPaths.foundDFC2 + "_top_ranked_hubs_pvalues", "F:\\\\Dropbox\\\\_Работа\\\\Programms\\\\DATA\\\\PhD\\\\Systems\\\\63mirnas\\\\output\\\\expression.txt", fn, FoldersPaths.foundDFC2 + "_top_ranked_hubs_pvalues_exp");
        //outils.add_aggregated_pvalues_to_pathways(all,FoldersPaths.foundDFC2 , "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\expression.txt", fn, FoldersPaths.foundDFC2 + "_aggr_pvalues",dbutils.hugo_by_id);
        //putils.get_most_confident_pathways(FoldersPaths.foundDFC2 + "_aggr_pvalues", FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked");

        // putils.add_connectivity_to_pathways(all, FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked", FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked_with_connectivity");
        //
        //outils.calculate_centrality_scores_for_paths(all, FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked", FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked__unique_overrepr", 1, 4);
        //outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked__unique_overrepr", FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        //outils.calculate_paths_degree(FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked_with_connectivity", FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked__unique_overrepr_w_names_pthconn");
        //putils.filterPathways_by_centrality_score(FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked__unique_overrepr_w_names_pthconn_flt", 3, 10000);
        //rand.calculate_p_values_phenotype_label_permutation(all, nutils, putils, outils, dbutils, FoldersPaths.C2random_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFC2 + "_aggr_pvalues_top_ranked__unique_overrepr_w_names_pthconn_flt", FoldersPaths.C2random_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        initial_hitlist_size = nutils.hg.size();
        rand.calculate_p_values_phenotype_label_permutation_nodes(all, nutils, putils, outils, dbutils, FoldersPaths.M63random_path + "foundDF.txt_top_ranked", FoldersPaths.foundDF63 + "_top_ranked_hubs", FoldersPaths.M63random_path + "RandomHitGenes_lists.txt", initial_hitlist_size, dbutils.hugo_by_id);

    }

    /**
     * Part 3 of pipeline for miRNA muscle differentiation screen in mixed
     * network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void mirna63Sys_part3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        putils.create_cyto_for_pathways(FoldersPaths.foundDF63 + "_top_ranked", FoldersPaths.outputF63 + "networks\\paths.txt", FoldersPaths.outputF63 + "networks\\", all, nutils.hg, nutils.fpl);

        nutils.load_hitlist_and_finalimpl(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        putils.filterPathways_by_length(FoldersPaths.foundDF63, "_max_4", 0, 4);
        //outils.find_overreprs(3, 2, FoldersPaths.foundDF63 + "_max_4", FoldersPaths.foundDF63 + "_max_4_overrepr");
        putils.create_cyto_for_paths(FoldersPaths.foundDF63 + "_max_4_overrepr", all, nutils.hg, nutils.fpl);
        putils.filterPathways(FoldersPaths.foundDF63, "_mir_in_the_middle", all, nutils.hg, nutils.fpl, "all", "", "");
        //outils.find_overreprs(4, 3, FoldersPaths.foundDF63 + "_mir_in_the_middle", FoldersPaths.foundDF63 + "_mir_in_the_middle_overrepr");
        putils.create_cyto_for_paths(FoldersPaths.foundDF63 + "_mir_in_the_middle_overrepr", all, nutils.hg, nutils.fpl);
        putils.find_miRNAs_on_pathways(FoldersPaths.foundDF63, FoldersPaths.foundDF63_mirnas, dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "");
        putils.find_miRNAs_on_pathways(FoldersPaths.foundDF63 + "_mir_in_the_middle", FoldersPaths.foundDF63_mirnas + "_mir_in_the_middle", dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "countonlymiddle");
        String[] mirnas = {"hsa-mir-125b", "hsa-mir-17", "hsa-let-7a"};
        for (String mirna : mirnas) {
            putils.filterPathways_2(FoldersPaths.foundDF63, "_" + mirna.substring(8), all, nutils.hg, nutils.fpl, mirna, 0, 4, 4, 6);
            //outils.find_overreprs(5, 3, FoldersPaths.foundDF63 + "_" + mirna.substring(8), FoldersPaths.foundDF63 + "_" + mirna.substring(8) + "_overrepr");
            putils.create_cyto_for_paths(FoldersPaths.foundDF63 + "_" + mirna.substring(8) + "_overrepr", all, nutils.hg, nutils.fpl);
        }
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        outils.filter_paths(FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names", "_min_occ_8", all, nutils.hg, nutils.fpl, 8, "occ");
        outils.filter_paths(FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names_min_occ_8", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.find_miRNAs_on_paths(FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names_min_occ_8_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 8, "");
        outils.filter_paths(FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.find_miRNAs_on_paths(FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 0, "");
        //
        outils.filter_paths(FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names", "_min_occ_20", all, nutils.hg, nutils.fpl, 20, "occ");
        outils.filter_paths(FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_min_occ_20", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.find_miRNAs_on_paths(FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_min_occ_20_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 20, "");
        outils.filter_paths(FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.find_miRNAs_on_paths(FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 0, "");
        putils.create_cyto_for_paths(FoldersPaths.foundDF63 + "_top_ranked_overrepr", all, nutils.hg, nutils.fpl);

    }

    

    /**
     * Part 1 of pipeline for transcriptome profiling muscle differentiation
     * screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void transcriptomic_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);
        nutils.find_pathway_for_list_BF_algorithm(all, dbutils.hugo, 6, FoldersPaths.foundDFC2, "2C");
        putils.find_the_shortest_paths(FoldersPaths.foundDFC2, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     * Part 2 of pipeline for transcriptome profiling muscle differentiation
     * screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void transcriptomic_part2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.add_connectivity_to_pathways(all, FoldersPaths.foundDFC2 + "_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked_with_connectivity");
        //
        outils.calculate_centrality_scores_for_paths(all, FoldersPaths.foundDFC2 + "_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr", 1, 4);
        outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.calculate_paths_degree(FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFC2 + "_top_ranked_with_connectivity", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names_pthconn");
    }

    /**
     * Part 3 of pipeline for transcriptome profiling muscle differentiation
     * screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void transcriptomic_part3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.filterPathways(FoldersPaths.foundDFC2 + "_not_flt", FoldersPaths.foundDFC2, all, nutils.hg, nutils.fpl, "mir_middle", "", "");
        putils.find_the_shortest_paths(FoldersPaths.foundDFC2, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
        putils.add_connectivity_to_pathways(all, FoldersPaths.foundDFC2 + "_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked_with_connectivity");
        //
        outils.calculate_centrality_scores_for_paths(all, FoldersPaths.foundDFC2 + "_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr", 1, 4);
        outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.calculate_paths_degree(FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFC2 + "_top_ranked_with_connectivity", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names_pthconn");

        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);
        putils.compare_two_paths_files(all, dbutils.hugo_by_id, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked__unique_overrepr_w_names_pthconn_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked__unique_overrepr_w_names_pthconn_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\c2_63.txt");
        putils.compare_two_nodes_files(all, dbutils.hugo_by_id, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked_hubs_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked_hubs_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\c2_63_hubs.txt");
        putils.find_hitgenes_on_paths(all, nutils, dbutils.hugo_by_id, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked__unique_overrepr_w_names_pthconn_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\c2_hits_63_overreps.txt");
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        putils.find_hitgenes_on_paths(all, nutils, dbutils.hugo_by_id, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked__unique_overrepr_w_names_pthconn_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\63_hits_c2_overreps.txt");

        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);

        outils.filter_paths(FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names", "_min_occ_10", all, nutils.hg, nutils.fpl, 10, "occ");
        outils.filter_paths(FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names_min_occ_10", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.find_miRNAs_on_paths(FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names_min_occ_10_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 10, "");
        outils.filter_paths(FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.find_miRNAs_on_paths(FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 0, "");
        //
        outils.filter_paths(FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr_w_names", "_min_occ_40", all, nutils.hg, nutils.fpl, 40, "occ");
        outils.filter_paths(FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_min_occ_40", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.find_miRNAs_on_paths(FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_min_occ_40_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 40, "");
        outils.filter_paths(FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr_w_names", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.find_miRNAs_on_paths(FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 0, "");
        //
        putils.create_cyto_for_paths(FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr", all, nutils.hg, nutils.fpl);
        putils.create_cyto_for_paths(FoldersPaths.foundDFC2 + "_top_ranked_overrepr", all, nutils.hg, nutils.fpl);
        putils.find_miRNAs_on_pathways(FoldersPaths.foundDFC2, FoldersPaths.foundDFC2_mirnas, dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "");

    }

    /**
     * Part 'p-values' of pipeline for transcriptome profiling muscle
     * differentiation screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void transcriptomic_random(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list = new ArrayList(nutils.hg.keySet());
        rand.permute_phenotype_label(all, nutils, putils, outils, dbutils, FoldersPaths.C2random_path, "C2R", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpC2);
        rand.calculate_p_values_phenotype_label_permutation(all, nutils, putils, outils, dbutils, FoldersPaths.C2random_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.C2random_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);
        rand.calculate_p_values_phenotype_label_permutation_nodes(all, nutils, putils, outils, dbutils, FoldersPaths.C2random_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked_hubs", FoldersPaths.C2random_path + "RandomHitGenes_lists.txt", initial_hitlist_size, dbutils.hugo_by_id);
    }

}
