package masterPATH;

import java.io.BufferedReader;
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

    /**
     * Add low confident connections for OGG1 protein
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     * @param hugo HGNC ids map
     */
    public Network ogg1_add_connections(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, Map<String, String[]> hugo) throws IOException {
        //hg.put(hugo.get(s)[0], s)
        dbutils.loadHIPPIE(false);
        String ogg1_id = hugo.get("OGG1")[0];
        String pogg1_id = "p" + hugo.get("OGG1")[0];
        //System.out.println(ogg1_id);        
        Network ogg1_con = new Network();
        ArrayList<Network> networks = new ArrayList();
        for (Network n : new Network[]{dbutils.hippie}) {
            for (String s : n.interactions.keySet()) {
                if (n.interactions.get(s).int1.id.equals(ogg1_id) || n.interactions.get(s).int2.id.equals(ogg1_id) || n.interactions.get(s).int1.id.equals(pogg1_id) || n.interactions.get(s).int2.id.equals(pogg1_id)) {
                    // ogg1_con = new Network();
                    ogg1_con.interactions.put(s, n.interactions.get(s));
                    System.out.println(s);
                    ogg1_con.nodes.put(n.interactions.get(s).int1.id, n.nodes.get(n.interactions.get(s).int1.id));
                    ogg1_con.nodes.put(n.interactions.get(s).int2.id, n.nodes.get(n.interactions.get(s).int2.id));
                }
            }
        }
        System.out.println("ogg1 network interactions " + ogg1_con.interactions.size() + "ogg1 network nodes " + ogg1_con.nodes.size());
        networks.add(ogg1_con);
        networks.add(all);
        Network all2 = nutils.merge_list_of_networks(networks);
        nutils.build_map_of_neighbors(all2);
        nutils.add_missing_genes_and_products(all2);
        nutils.build_map_of_neighbors(all2);
        return all2;
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
     * Part 'pvalues' of pipeline for DNA repair screen in ppi network
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
     * Part 'random network' of pipeline for DNA repair screen in ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void ogg1_random2_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlOG, FoldersPaths.fpOG, dbutils.hugo);

        rand.calculate_p_values_nodes_random_networks_ppi(all, nutils, putils, outils, dbutils, FoldersPaths.OGrandom_path2_ppi + "foundDF.txt", 1000, FoldersPaths.foundDFOG_ppi + "_top_ranked_hubs", "");
    }

    /**
     * Part 1 of pipeline for DNA repair screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void ogg1_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlOG, FoldersPaths.fpOG, dbutils.hugo);
        nutils.find_pathway_for_list_BF_algorithm(all, dbutils.hugo, 5, FoldersPaths.foundDFOG, "OG");
        putils.find_the_shortest_paths(FoldersPaths.foundDFOG, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     * Part 2 of pipeline for DNA repair screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void ogg1_part2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.add_connectivity_to_pathways(all, FoldersPaths.foundDFOG + "_top_ranked", FoldersPaths.foundDFOG + "_top_ranked_with_connectivity");
        //
        outils.calculate_centrality_scores_for_paths(all, FoldersPaths.foundDFOG + "_top_ranked", FoldersPaths.foundDFOG + "_top_ranked__unique_overrepr", 2, 3);
        outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFOG + "_top_ranked__unique_overrepr", FoldersPaths.foundDFOG + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.calculate_paths_degree(FoldersPaths.foundDFOG + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFOG + "_top_ranked_with_connectivity", FoldersPaths.foundDFOG + "_top_ranked__unique_overrepr_w_names_pthconn");

        putils.filter_pathways_by_connectivity(FoldersPaths.foundDFOG + "_top_ranked_with_connectivity", FoldersPaths.foundDFOG + "_top_ranked_flt_con", FoldersPaths.foundDFOG + "_top_ranked_flt_con_for_overreps", 0, 60, 0, 60);
        outils.calculate_centrality_scores_for_paths(all, FoldersPaths.foundDFOG + "_top_ranked_flt_con_for_overreps", FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr", 2, 4);
        outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr", FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.calculate_paths_degree(FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_w_names", FoldersPaths.foundDFOG + "_top_ranked_with_connectivity", FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_w_names_pthconn");
    }

    /**
     * Part 3 of pipeline for DNA repair screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void ogg1_part3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        for (int i = 0; i <= 1200; i = i + 100) {
            putils.filter_pathways_by_connectivity(FoldersPaths.foundDFOG + "_top_ranked_with_connectivity", FoldersPaths.foundDFOG + "_top_ranked_flt_con_" + i, FoldersPaths.foundDFOG + "_top_ranked_flt_con_for_overreps_" + i, 0, i, 0, i);
            outils.calculate_centrality_scores_for_paths(all, FoldersPaths.foundDFOG + "_top_ranked_flt_con_for_overreps_" + i, FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_" + i, 2, 4);
            outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_" + i, FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_w_names_" + i, all, dbutils.hugo_by_id);
            outils.calculate_paths_degree(FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_w_names_" + i, FoldersPaths.foundDFOG + "_top_ranked_with_connectivity", FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_w_names_pthconn_" + i);
        }
    }

    /**
     * Part 4 of pipeline for DNA repair screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void ogg1_part4(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlOG, FoldersPaths.fpOG, dbutils.hugo);
        putils.create_cyto_for_pathways(FoldersPaths.foundDFOG + "_top_ranked", FoldersPaths.outputFOG + "networks\\paths.txt", FoldersPaths.outputFOG + "networks\\", all, nutils.hg, nutils.fpl);
    }

    /**
     * Part 'p-values' of pipeline for DNA repair screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void ogg1_random(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlOG, FoldersPaths.fpOG, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list = new ArrayList(nutils.hg.keySet());
        rand.permute_phenotype_label(all, nutils, putils, outils, dbutils, FoldersPaths.OGrandom_path, "OGR", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpOG);
        rand.calculate_p_values_phenotype_label_permutation(all, nutils, putils, outils, dbutils, FoldersPaths.OGrandom_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFOG + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.OGrandom_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 5, dbutils.hugo_by_id);
    }

    /**
     * Part 1 of pipeline for arp2/3 screen in ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void arp_part1_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlARP, FoldersPaths.fpARP_ppi, dbutils.hugo);
        nutils.find_pathway_for_list_BF_algorithm_ppi(all, dbutils.hugo, 5, FoldersPaths.foundDFARP_ppi, "ARP");
        putils.find_the_shortest_paths(FoldersPaths.foundDFARP_ppi, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     * Part 2 of pipeline for arp2/3 screen in ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void arp_part2_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.add_connectivity_to_pathways(all, FoldersPaths.foundDFARP_ppi + "_top_ranked", FoldersPaths.foundDFARP_ppi + "_top_ranked_with_connectivity");
        //
        outils.calculate_centrality_scores_for_paths_ppi(all, FoldersPaths.foundDFARP_ppi + "_top_ranked", FoldersPaths.foundDFARP_ppi + "_top_ranked__unique_overrepr", 2, 3);
        outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFARP_ppi + "_top_ranked__unique_overrepr", FoldersPaths.foundDFARP_ppi + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.calculate_paths_degree(FoldersPaths.foundDFARP_ppi + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFARP_ppi + "_top_ranked_with_connectivity", FoldersPaths.foundDFARP_ppi + "_top_ranked__unique_overrepr_w_names_pthconn");
    }

    /**
     * Part 'p-values' of pipeline for arp2/3 screen in ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void arp_random_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlARP, FoldersPaths.fpARP_ppi, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list = new ArrayList(nutils.hg.keySet());
        rand.permute_phenotype_label(all, nutils, putils, outils, dbutils, FoldersPaths.ARPrandom_path_ppi, "ARPR", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpARP_ppi);
        rand.calculate_p_values_phenotype_label_permutation_ppi(all, nutils, putils, outils, dbutils, FoldersPaths.ARPrandom_path_ppi + "foundDF.txt_top_ranked", FoldersPaths.foundDFARP_ppi + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.ARPrandom_path_ppi + "RandomHitGenes_lists.txt", initial_hitlist_size, 5, dbutils.hugo_by_id);
    }

    /**
     * Add low confident interactions for LIMCH1 protein
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public Network arp_add_connections(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, Map<String, String[]> hugo, String hitlist) throws IOException {
        String t, limch1_id, plimch1_id;
        Network limch1_con = new Network();
        ArrayList<Network> networks = new ArrayList();

        BufferedReader rd = new BufferedReader(new FileReader(hitlist));
        dbutils.loadHIPPIE(false);
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlARP, FoldersPaths.fpARP, dbutils.hugo);
        while ((t = rd.readLine()) != null) {
            limch1_id = hugo.get(t)[0];
            plimch1_id = "p" + hugo.get(t)[0];
            for (Network n : new Network[]{dbutils.hippie}) {
                for (String s : n.interactions.keySet()) {
                    if (n.interactions.get(s).int1.id.equals(limch1_id) || n.interactions.get(s).int2.id.equals(limch1_id) || n.interactions.get(s).int1.id.equals(plimch1_id) || n.interactions.get(s).int2.id.equals(plimch1_id)) {
                        limch1_con.interactions.put(s, n.interactions.get(s));
                        //System.out.println(s);
                        limch1_con.nodes.put(n.interactions.get(s).int1.id, n.nodes.get(n.interactions.get(s).int1.id));
                        limch1_con.nodes.put(n.interactions.get(s).int2.id, n.nodes.get(n.interactions.get(s).int2.id));
                    }
                }
            }
        }

        System.out.println("limch1 network interactions " + limch1_con.interactions.size() + " limch1 network nodes " + limch1_con.nodes.size());
        networks.add(limch1_con);
        networks.add(all);
        Network all2 = nutils.merge_list_of_networks(networks);
        nutils.build_map_of_neighbors(all2);
        nutils.add_missing_genes_and_products(all2);
        nutils.build_map_of_neighbors(all2);
        return all2;
    }

    /**
     * Part 1 of pipeline for arp2/3 screen in mixed network Hit list only LMCH1
     * gene
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void arp_part1_LIMCH1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlARP_LIMCH1, FoldersPaths.fpARP, dbutils.hugo);
        nutils.find_pathway_for_list_BF_algorithm(all, dbutils.hugo, 5, FoldersPaths.foundDFARP + "LIMCH1", "ARPLIMCH1");
        putils.find_the_shortest_paths(FoldersPaths.foundDFARP + "LIMCH1", "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);

        putils.create_cyto_for_pathways(FoldersPaths.foundDFARP + "LIMCH1_top_ranked", FoldersPaths.outputFARP + "networks\\paths.txt", FoldersPaths.outputFARP + "networks\\", all, nutils.hg, nutils.fpl);
    }

    /**
     * Part 1 of pipeline for arp2/3 screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void arp_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlARP, FoldersPaths.fpARP, dbutils.hugo);
        nutils.find_pathway_for_list_BF_algorithm(all, dbutils.hugo, 5, FoldersPaths.foundDFARP, "ARP");
        putils.find_the_shortest_paths(FoldersPaths.foundDFARP, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     * Part 2 of pipeline for arp2/3 screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void arp_part2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.add_connectivity_to_pathways(all, FoldersPaths.foundDFARP + "_top_ranked", FoldersPaths.foundDFARP + "_top_ranked_with_connectivity");

        outils.calculate_centrality_scores_for_paths(all, FoldersPaths.foundDFARP + "_top_ranked", FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr", 2, 4);
        outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr", FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.calculate_paths_degree(FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFARP + "_top_ranked_with_connectivity", FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr_w_names_pthconn");

        putils.filter_pathways_by_connectivity(FoldersPaths.foundDFARP + "_top_ranked_with_connectivity", FoldersPaths.foundDFARP + "_top_ranked_flt_con", FoldersPaths.foundDFARP + "_top_ranked_flt_con_for_overreps", 0, 60, 0, 60);
        outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFARP + "_top_ranked_flt_con_unique_overrepr", FoldersPaths.foundDFARP + "_top_ranked_flt_con_unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.calculate_paths_degree(FoldersPaths.foundDFARP + "_top_ranked_flt_con_unique_overrepr_w_names", FoldersPaths.foundDFARP + "_top_ranked_with_connectivity", FoldersPaths.foundDFARP + "_top_ranked_flt_con_unique_overrepr_w_names_pthconn");

    }

    /**
     * Part 'p-values' of pipeline for arp2/3 screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void arp_random(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlARP, FoldersPaths.fpARP, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list = new ArrayList(nutils.hg.keySet());
        rand.permute_phenotype_label(all, nutils, putils, outils, dbutils, FoldersPaths.ARPrandom_path, "APR", 5000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpARP);
        rand.calculate_p_values_phenotype_label_permutation(all, nutils, putils, outils, dbutils, FoldersPaths.ARPrandom_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.ARPrandom_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);

        rand.calculate_p_values_paths_random_networks(all, nutils, putils, outils, dbutils, FoldersPaths.ARPrandom_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.ARPrandom_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);
    }

    /**
     * Part 3 of pipeline for arp2/3 screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void arp_part3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlARP, FoldersPaths.fpARP, dbutils.hugo);
        putils.create_cyto_for_pathways(FoldersPaths.foundDFARP + "_top_ranked", FoldersPaths.outputFARP + "paths.txt", FoldersPaths.outputFARP + "", all, nutils.hg, nutils.fpl);
    }

    /**
     * Filter pathways for arp2/3 screen in mixed network
     *
     * @param all Network
     */
    public void arp_filter(Network all) throws FileNotFoundException, IOException {
        String s, file = "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\ARP\\input\\genes_filter_out.txt";
        String[] ss;
        int count = 0, count_top = 0;
        Map<String, String> genes_filter = new HashMap();
        BufferedReader rd = new BufferedReader(new FileReader(file));
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            genes_filter.put(ss[1], ss[0]);
        }
        BufferedReader rd_pathways = new BufferedReader(new FileReader(FoldersPaths.foundDFARP));
        while ((s = rd_pathways.readLine()) != null) {
            ss = s.split("\t");
            for (int i = 2; i < ss.length - 1; i++) {
                if (genes_filter.containsKey(all.interactions.get(ss[i]).int1.id) || genes_filter.containsKey(all.interactions.get(ss[i]).int2.id)) {
                    count++;
                }
            }
        }
        BufferedReader rd_top_pathways = new BufferedReader(new FileReader(FoldersPaths.foundDFARP + "_top_ranked"));
        while ((s = rd_top_pathways.readLine()) != null) {
            ss = s.split("\t");
            for (int i = 2; i < ss.length - 1; i++) {
                if (genes_filter.containsKey(all.interactions.get(ss[i]).int1.id) || genes_filter.containsKey(all.interactions.get(ss[i]).int2.id)) {
                    count_top++;
                }
            }
        }
        rd.close();
        rd_pathways.close();
        rd_top_pathways.close();
        System.out.println(count);
        System.out.println(count_top);
    }

    /**
     * Part 1 of pipeline for adenocarcinoma screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void adenocarcinoma_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlAd, FoldersPaths.fpAd, dbutils.hugo);
        nutils.find_pathway_for_list_BF_algorithm(all, dbutils.hugo, 5, FoldersPaths.foundDFAd, "AD");
        putils.find_the_shortest_paths(FoldersPaths.foundDFAd, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     * Part 2 of pipeline for adenocarcinoma screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void adenocarcinoma_part2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.add_connectivity_to_pathways(all, FoldersPaths.foundDFAd + "_top_ranked", FoldersPaths.foundDFAd + "_top_ranked_with_connectivity");
        //
        outils.calculate_centrality_scores_for_paths(all, FoldersPaths.foundDFAd + "_top_ranked", FoldersPaths.foundDFAd + "_top_ranked__unique_overrepr", 2, 4);
        outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFAd + "_top_ranked__unique_overrepr", FoldersPaths.foundDFAd + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.calculate_paths_degree(FoldersPaths.foundDFAd + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFAd + "_top_ranked_with_connectivity", FoldersPaths.foundDFAd + "_top_ranked__unique_overrepr_w_names_pthconn");
    }

    /**
     * Part 'p-values' of pipeline for adenocarcinoma screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void adenocarcinoma_random(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlAd, FoldersPaths.fpAd, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list = new ArrayList(nutils.hg.keySet());
        rand.permute_phenotype_label(all, nutils, putils, outils, dbutils, FoldersPaths.Adrandom_path, "AdR", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpAd);
        rand.calculate_p_values_phenotype_label_permutation(all, nutils, putils, outils, dbutils, FoldersPaths.Adrandom_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFAd + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.Adrandom_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);
    }

    /**
     * Part 1 of pipeline for adenocarcinoma screen in ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void adenocarcinoma_part1_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlAd, FoldersPaths.fpAd_ppi, dbutils.hugo);
        nutils.find_pathway_for_list_BF_algorithm_ppi(all, dbutils.hugo, 5, FoldersPaths.foundDFAd_ppi, "Ad");
        putils.find_the_shortest_paths(FoldersPaths.foundDFAd_ppi, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     * Part 2 of pipeline for adenocarcinoma screen in ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void adenocarcinoma_part2_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.add_connectivity_to_pathways(all, FoldersPaths.foundDFAd_ppi + "_top_ranked", FoldersPaths.foundDFAd_ppi + "_top_ranked_with_connectivity");
        outils.calculate_centrality_scores_for_paths_ppi(all, FoldersPaths.foundDFAd_ppi + "_top_ranked", FoldersPaths.foundDFAd_ppi + "_top_ranked__unique_overrepr", 2, 4);
        outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFAd_ppi + "_top_ranked__unique_overrepr", FoldersPaths.foundDFAd_ppi + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.calculate_paths_degree(FoldersPaths.foundDFAd_ppi + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFAd_ppi + "_top_ranked_with_connectivity", FoldersPaths.foundDFAd_ppi + "_top_ranked__unique_overrepr_w_names_pthconn");
    }

    /**
     * Part 'p-values' of pipeline for adenocarcinoma screen in ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void adenocarcinoma_random_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlAd, FoldersPaths.fpAd_ppi, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list = new ArrayList(nutils.hg.keySet());
        rand.permute_phenotype_label(all, nutils, putils, outils, dbutils, FoldersPaths.Adrandom_path_ppi, "AdR", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpAd_ppi);
        rand.calculate_p_values_phenotype_label_permutation_ppi(all, nutils, putils, outils, dbutils, FoldersPaths.Adrandom_path_ppi + "foundDF.txt_top_ranked", FoldersPaths.foundDFAd_ppi + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.Adrandom_path_ppi + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);
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
    }

    /**
     * Part 3 + filter of pipeline for miRNA muscle differentiation screen in
     * mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public static void mirna63Sys_part3_3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
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
    }

    /**
     * Part 4 of pipeline for miRNA muscle differentiation screen in mixed
     * network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public static void mirna63Sys_part4(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
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
    public void case2_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
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
    public void case2_part2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
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
    public void case2_part3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.filterPathways(FoldersPaths.foundDFC2 + "_not_flt", FoldersPaths.foundDFC2, all, nutils.hg, nutils.fpl, "mir_middle", "", "");
        putils.find_the_shortest_paths(FoldersPaths.foundDFC2, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
        putils.add_connectivity_to_pathways(all, FoldersPaths.foundDFC2 + "_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked_with_connectivity");
        //
        outils.calculate_centrality_scores_for_paths(all, FoldersPaths.foundDFC2 + "_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr", 1, 4);
        outils.add_hgnc_symbols_to_paths(FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.calculate_paths_degree(FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFC2 + "_top_ranked_with_connectivity", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names_pthconn");
    }

    /**
     * Part 4 of pipeline for transcriptome profiling muscle differentiation
     * screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void case2_part4(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);
        putils.compare_two_paths_files(all, dbutils.hugo_by_id, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked__unique_overrepr_w_names_pthconn_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked__unique_overrepr_w_names_pthconn_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\c2_63.txt");
        putils.compare_two_nodes_files(all, dbutils.hugo_by_id, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked_hubs_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked_hubs_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\c2_63_hubs.txt");
        putils.find_hitgenes_on_paths(all, nutils, dbutils.hugo_by_id, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked__unique_overrepr_w_names_pthconn_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\c2_hits_63_overreps.txt");
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        putils.find_hitgenes_on_paths(all, nutils, dbutils.hugo_by_id, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked__unique_overrepr_w_names_pthconn_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\63_hits_c2_overreps.txt");
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
    public void case2_random(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list = new ArrayList(nutils.hg.keySet());
        rand.permute_phenotype_label(all, nutils, putils, outils, dbutils, FoldersPaths.C2random_path, "C2R", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpC2);
        rand.calculate_p_values_phenotype_label_permutation(all, nutils, putils, outils, dbutils, FoldersPaths.C2random_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.C2random_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);
        rand.calculate_p_values_phenotype_label_permutation_nodes(all, nutils, putils, outils, dbutils, FoldersPaths.C2random_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked_hubs", FoldersPaths.C2random_path + "RandomHitGenes_lists.txt", initial_hitlist_size, dbutils.hugo_by_id);
    }

    /**
     * Part 'p-values on random networks' of pipeline for transcriptome
     * profiling muscle differentiation screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void case2_random2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.load_hitlist_and_finalimpl(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);
        rand.calculate_p_values_nodes_random_networks(all, nutils, putils, outils, dbutils, FoldersPaths.C2random_path2 + "foundDF.txt", 1, FoldersPaths.foundDFC2 + "_top_ranked_hubs", "");
    }

    /**
     * Part 2 + filter of pipeline for transcriptome profiling muscle
     * differentiation screen in mixed network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public static void case2_part2_2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
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

}

/*
 public static void muscleDifSys_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
 nutils.load_hitlist_and_finalimpl(FoldersPaths.hl, FoldersPaths.fp, dbutils.hugo);
 //      nutils.getLongPathsforListDF(all, dbutils.hugo, 6, FoldersPaths.foundDF, "HP");
 putils.find_the_shortest_paths(FoldersPaths.foundDF, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
 putils.find_the_shortest_paths(FoldersPaths.foundDF, "_top_ranked_ppi_min_plus_1", all, nutils.hg, nutils.fpl, 1, 1, 0, 0);
 outils.find_overreprs(4, 2, FoldersPaths.foundDF, FoldersPaths.foundDF + "_overrepr");
 outils.find_overreprs(4, 2, FoldersPaths.foundDF + "_top_ranked", FoldersPaths.foundDF + "_top_ranked_overrepr");
 outils.find_overreprs(4, 2, FoldersPaths.foundDF + "_top_ranked_ppi_min_plus_1", FoldersPaths.foundDF + "_top_ranked_ppi_min_plus_1_overrepr");
 outils.put_gene_symbols_for_paths(FoldersPaths.foundDF + "_overrepr", FoldersPaths.foundDF + "_overrepr_w_names", all, dbutils.hugo_by_id);
 outils.put_gene_symbols_for_paths(FoldersPaths.foundDF + "_top_ranked_overrepr", FoldersPaths.foundDF + "_top_ranked_overrepr_w_names", all, dbutils.hugo_by_id);
 outils.put_gene_symbols_for_paths(FoldersPaths.foundDF + "_top_ranked_ppi_min_plus_1_overrepr", FoldersPaths.foundDF + "_top_ranked_ppi_min_plus_1_overrepr_w_names", all, dbutils.hugo_by_id);
 }

 public static void muscleDifSys_part3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
 nutils.load_hitlist_and_finalimpl(FoldersPaths.hl, FoldersPaths.fp, dbutils.hugo);
 putils.find_miRNAs_on_pathways(FoldersPaths.foundDF, FoldersPaths.foundDF_mirnas, dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "");
 putils.find_miRNAs_on_pathways(FoldersPaths.foundDF, FoldersPaths.foundDF_mirnas, dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "countonlymiddle");
 putils.find_miRNAs_on_pathways(FoldersPaths.foundDF + "_top_ranked", FoldersPaths.foundDF_mirnas + "_top_ranked", dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "");
 putils.find_miRNAs_on_pathways(FoldersPaths.foundDF + "_top_ranked_ppi_min_plus_1", FoldersPaths.foundDF_mirnas + "_top_ranked_ppi_min_plus_1", dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "");
 String[] mirnas = {"hsa-mir-34a", "hsa-mir-200b", "hsa-mir-27a", "hsa-mir-449a", "hsa-mir-30a", "hsa-mir-148a", "hsa-mir-145", "hsa-mir-206", "hsa-mir-125b", "hsa-mir-125b", "hsa-let-7a"};
 for (String mirna : mirnas) {
 putils.filterPathways(FoldersPaths.foundDF, "_" + mirna.substring(8), all, nutils.hg, nutils.fpl, "mirna", mirna, "miRNA");
 outils.find_overreprs(5, 3, FoldersPaths.foundDF + "_" + mirna.substring(8), FoldersPaths.foundDF + "_" + mirna.substring(8) + "_overrepr");
 putils.create_cyto_for_paths(FoldersPaths.foundDF + "_" + mirna.substring(8) + "_overrepr", all, nutils.hg, nutils.fpl);
 }
 putils.filterPathways_by_length(FoldersPaths.outputF + "foundDF_28_11.txt_206", "_max_5", 0, 5);
 outils.find_overreprs(4, 3, FoldersPaths.outputF + "foundDF_28_11.txt_206_max_5", FoldersPaths.outputF + "foundDF_28_11.txt_206_max_5_overrepr");
 putils.create_cyto_for_paths(FoldersPaths.outputF + "foundDF_28_11.txt_206_max_5_overrepr", all, nutils.hg, nutils.fpl);
 putils.filterPathways(FoldersPaths.outputF + "foundDF_28_11.txt_206", "_filter_mirna", all, nutils.hg, nutils.fpl, "all", "", "");
 outils.find_overreprs(4, 3, FoldersPaths.outputF + "foundDF_28_11.txt_206_filter_mirna", FoldersPaths.outputF + "foundDF_28_11.txt_206_filter_mirna_overrepr");
 putils.create_cyto_for_paths(FoldersPaths.outputF + "foundDF_28_11.txt_206_filter_mirna_overrepr", all, nutils.hg, nutils.fpl);
 }

 public static void muscleDifSys_part4(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
 nutils.load_hitlist_and_finalimpl(FoldersPaths.hl, FoldersPaths.fp, dbutils.hugo);
 putils.find_miRNAs_on_paths(FoldersPaths.foundDF + "_overrepr_w_names", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 3, "");
 putils.find_miRNAs_on_paths(FoldersPaths.foundDF + "_top_ranked_overrepr_w_names", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 3, "");
 putils.find_miRNAs_on_paths(FoldersPaths.foundDF + "_top_ranked_ppi_min_plus_1_overrepr_w_names", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 3, "");
 }

 */
