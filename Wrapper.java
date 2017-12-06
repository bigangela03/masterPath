/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package masterPATH;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author Natalia Rubanova
 */
public class Wrapper {

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param hugo
     * @return
     * @throws IOException
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
        nutils.buildNeighbours(all2);
        nutils.add_genes_and_products(all2);
        nutils.buildNeighbours(all2);
        return all2;
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void ogg1_part1_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlOG, FoldersPaths.fpOG_ppi, dbutils.hugo);
        nutils.getLongPathsforListBF_ppi(all, dbutils.hugo, 5, FoldersPaths.foundDFOG_ppi, "OG");
        putils.rank_pathways(FoldersPaths.foundDFOG_ppi, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void ogg1_part2_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.calculate_connectivity_for_pathways(all, FoldersPaths.foundDFOG_ppi + "_top_ranked", FoldersPaths.foundDFOG_ppi + "_top_ranked_with_connectivity");
        //
        outils.get_centrality_paths_wo_redundancy_ppi(all, FoldersPaths.foundDFOG_ppi + "_top_ranked", FoldersPaths.foundDFOG_ppi + "_top_ranked__unique_overrepr", 2, 3);
        outils.put_names_on_common_parts_wo_redundancy(FoldersPaths.foundDFOG_ppi + "_top_ranked__unique_overrepr", FoldersPaths.foundDFOG_ppi + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.put_connectivity_on_paths_wo_redundancy(FoldersPaths.foundDFOG_ppi + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFOG_ppi + "_top_ranked_with_connectivity", FoldersPaths.foundDFOG_ppi + "_top_ranked__unique_overrepr_w_names_pthconn");
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param rand
     * @throws IOException
     */
    public void ogg1_random_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlOG, FoldersPaths.fpOG_ppi, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list=new ArrayList(nutils.hg.keySet());
        rand.build_paths_for_random_hit_lists(all, nutils, putils, outils, dbutils, FoldersPaths.OGrandom_path_ppi, "OGR", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpOG_ppi);
        rand.calculate_p_values_ppi(all, nutils, putils, outils, dbutils, FoldersPaths.OGrandom_path_ppi + "foundDF.txt_top_ranked", FoldersPaths.foundDFOG_ppi + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.OGrandom_path_ppi + "RandomHitGenes_lists.txt", initial_hitlist_size, 5, dbutils.hugo_by_id);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void ogg1_part4_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlOG, FoldersPaths.fpOG_ppi, dbutils.hugo);
        putils.create_cyto_for_list(FoldersPaths.foundDFOG_ppi + "_top_ranked", FoldersPaths.outputFOG + "ppi\\networks\\paths.txt", FoldersPaths.outputFOG + "ppi\\networks\\", all, nutils.hg, nutils.fpl);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param rand
     * @throws IOException
     */
    public void ogg1_random2_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlOG, FoldersPaths.fpOG, dbutils.hugo);

        rand.test_on_random_networks_ppi(all, nutils, putils, outils, dbutils, FoldersPaths.OGrandom_path2_ppi + "foundDF.txt", 1000, FoldersPaths.foundDFOG_ppi + "_top_ranked_hubs", "");
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void ogg1_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlOG, FoldersPaths.fpOG, dbutils.hugo);
        nutils.getLongPathsforListBF(all, dbutils.hugo, 5, FoldersPaths.foundDFOG, "OG");
        putils.rank_pathways(FoldersPaths.foundDFOG, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void ogg1_part2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.calculate_connectivity_for_pathways(all, FoldersPaths.foundDFOG + "_top_ranked", FoldersPaths.foundDFOG + "_top_ranked_with_connectivity");
        //
        outils.get_centrality_paths_wo_redundancy(all, FoldersPaths.foundDFOG + "_top_ranked", FoldersPaths.foundDFOG + "_top_ranked__unique_overrepr", 2, 3);
        outils.put_names_on_common_parts_wo_redundancy(FoldersPaths.foundDFOG + "_top_ranked__unique_overrepr", FoldersPaths.foundDFOG + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.put_connectivity_on_paths_wo_redundancy(FoldersPaths.foundDFOG + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFOG + "_top_ranked_with_connectivity", FoldersPaths.foundDFOG + "_top_ranked__unique_overrepr_w_names_pthconn");

        putils.filter_pathways_by_connectivity(FoldersPaths.foundDFOG + "_top_ranked_with_connectivity", FoldersPaths.foundDFOG + "_top_ranked_flt_con", FoldersPaths.foundDFOG + "_top_ranked_flt_con_for_overreps", 0, 60, 0, 60);
        outils.get_centrality_paths_wo_redundancy(all, FoldersPaths.foundDFOG + "_top_ranked_flt_con_for_overreps", FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr", 2, 4);
        outils.put_names_on_common_parts_wo_redundancy(FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr", FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.put_connectivity_on_paths_wo_redundancy(FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_w_names", FoldersPaths.foundDFOG + "_top_ranked_with_connectivity", FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_w_names_pthconn");
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void ogg1_part3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        for (int i = 0; i <= 1200; i = i + 100) {
            putils.filter_pathways_by_connectivity(FoldersPaths.foundDFOG + "_top_ranked_with_connectivity", FoldersPaths.foundDFOG + "_top_ranked_flt_con_" + i, FoldersPaths.foundDFOG + "_top_ranked_flt_con_for_overreps_" + i, 0, i, 0, i);
            outils.get_centrality_paths_wo_redundancy(all, FoldersPaths.foundDFOG + "_top_ranked_flt_con_for_overreps_" + i, FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_" + i, 2, 4);
            outils.put_names_on_common_parts_wo_redundancy(FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_" + i, FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_w_names_" + i, all, dbutils.hugo_by_id);
            outils.put_connectivity_on_paths_wo_redundancy(FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_w_names_" + i, FoldersPaths.foundDFOG + "_top_ranked_with_connectivity", FoldersPaths.foundDFOG + "_top_ranked_flt_con_unique_overrepr_w_names_pthconn_" + i);
        }
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void ogg1_part4(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlOG, FoldersPaths.fpOG, dbutils.hugo);
        putils.create_cyto_for_list(FoldersPaths.foundDFOG + "_top_ranked", FoldersPaths.outputFOG + "networks\\paths.txt", FoldersPaths.outputFOG + "networks\\", all, nutils.hg, nutils.fpl);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param rand
     * @throws IOException
     */
    public void ogg1_random(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlOG, FoldersPaths.fpOG, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list=new ArrayList(nutils.hg.keySet());
        rand.build_paths_for_random_hit_lists(all, nutils, putils, outils, dbutils, FoldersPaths.OGrandom_path, "OGR", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpOG);
        rand.calculate_p_values(all, nutils, putils, outils, dbutils, FoldersPaths.OGrandom_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFOG + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.OGrandom_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 5, dbutils.hugo_by_id);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void arp_part1_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlARP, FoldersPaths.fpARP_ppi, dbutils.hugo);
        nutils.getLongPathsforListBF_ppi(all, dbutils.hugo, 5, FoldersPaths.foundDFARP_ppi, "ARP");
        putils.rank_pathways(FoldersPaths.foundDFARP_ppi, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void arp_part2_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.calculate_connectivity_for_pathways(all, FoldersPaths.foundDFARP_ppi + "_top_ranked", FoldersPaths.foundDFARP_ppi + "_top_ranked_with_connectivity");
        //
        outils.get_centrality_paths_wo_redundancy_ppi(all, FoldersPaths.foundDFARP_ppi + "_top_ranked", FoldersPaths.foundDFARP_ppi + "_top_ranked__unique_overrepr", 2, 3);
        outils.put_names_on_common_parts_wo_redundancy(FoldersPaths.foundDFARP_ppi + "_top_ranked__unique_overrepr", FoldersPaths.foundDFARP_ppi + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.put_connectivity_on_paths_wo_redundancy(FoldersPaths.foundDFARP_ppi + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFARP_ppi + "_top_ranked_with_connectivity", FoldersPaths.foundDFARP_ppi + "_top_ranked__unique_overrepr_w_names_pthconn");
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param rand
     * @throws IOException
     */
    public void arp_random_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlARP, FoldersPaths.fpARP_ppi, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list=new ArrayList(nutils.hg.keySet());
        rand.build_paths_for_random_hit_lists(all, nutils, putils, outils, dbutils, FoldersPaths.ARPrandom_path_ppi, "ARPR", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpARP_ppi);
        rand.calculate_p_values_ppi(all, nutils, putils, outils, dbutils, FoldersPaths.ARPrandom_path_ppi + "foundDF.txt_top_ranked", FoldersPaths.foundDFARP_ppi + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.ARPrandom_path_ppi + "RandomHitGenes_lists.txt", initial_hitlist_size, 5, dbutils.hugo_by_id);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param hugo
     * @param hitlist
     * @return
     * @throws IOException
     */
    public Network arp_add_connections(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, Map<String, String[]> hugo, String hitlist) throws IOException {
        String t, limch1_id, plimch1_id;
        Network limch1_con = new Network();
        ArrayList<Network> networks = new ArrayList();

        BufferedReader rd = new BufferedReader(new FileReader(hitlist));
        dbutils.loadHIPPIE(false);
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlARP, FoldersPaths.fpARP, dbutils.hugo);
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
        nutils.buildNeighbours(all2);
        nutils.add_genes_and_products(all2);
        nutils.buildNeighbours(all2);
        return all2;
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void arp_part1_LIMCH1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlARP_LIMCH1, FoldersPaths.fpARP, dbutils.hugo);
        nutils.getLongPathsforListBF(all, dbutils.hugo, 5, FoldersPaths.foundDFARP + "LIMCH1", "ARPLIMCH1");
        putils.rank_pathways(FoldersPaths.foundDFARP + "LIMCH1", "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);

        putils.create_cyto_for_list(FoldersPaths.foundDFARP + "LIMCH1_top_ranked", FoldersPaths.outputFARP + "networks\\paths.txt", FoldersPaths.outputFARP + "networks\\", all, nutils.hg, nutils.fpl);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void arp_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlARP, FoldersPaths.fpARP, dbutils.hugo);
        nutils.getLongPathsforListBF(all, dbutils.hugo, 5, FoldersPaths.foundDFARP, "ARP");
        putils.rank_pathways(FoldersPaths.foundDFARP, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void arp_part2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.calculate_connectivity_for_pathways(all, FoldersPaths.foundDFARP + "_top_ranked", FoldersPaths.foundDFARP + "_top_ranked_with_connectivity");

        outils.get_centrality_paths_wo_redundancy(all, FoldersPaths.foundDFARP + "_top_ranked", FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr", 2, 4);
        outils.put_names_on_common_parts_wo_redundancy(FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr", FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.put_connectivity_on_paths_wo_redundancy(FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFARP + "_top_ranked_with_connectivity", FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr_w_names_pthconn");

        putils.filter_pathways_by_connectivity(FoldersPaths.foundDFARP + "_top_ranked_with_connectivity", FoldersPaths.foundDFARP + "_top_ranked_flt_con", FoldersPaths.foundDFARP + "_top_ranked_flt_con_for_overreps", 0, 60, 0, 60);
        outils.put_names_on_common_parts_wo_redundancy(FoldersPaths.foundDFARP + "_top_ranked_flt_con_unique_overrepr", FoldersPaths.foundDFARP + "_top_ranked_flt_con_unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.put_connectivity_on_paths_wo_redundancy(FoldersPaths.foundDFARP + "_top_ranked_flt_con_unique_overrepr_w_names", FoldersPaths.foundDFARP + "_top_ranked_with_connectivity", FoldersPaths.foundDFARP + "_top_ranked_flt_con_unique_overrepr_w_names_pthconn");
         
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param rand
     * @throws IOException
     */
    public void arp_random(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlARP, FoldersPaths.fpARP, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list=new ArrayList(nutils.hg.keySet());
        rand.build_paths_for_random_hit_lists(all, nutils, putils, outils, dbutils, FoldersPaths.ARPrandom_path, "APR", 5000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpARP);
        rand.calculate_p_values(all, nutils, putils, outils, dbutils, FoldersPaths.ARPrandom_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.ARPrandom_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);

        rand.calculate_p_values_random_networks(all, nutils, putils, outils, dbutils, FoldersPaths.ARPrandom_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFARP + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.ARPrandom_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void arp_part3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlARP, FoldersPaths.fpARP, dbutils.hugo);
        putils.create_cyto_for_list(FoldersPaths.foundDFARP + "_top_ranked", FoldersPaths.outputFARP + "paths.txt", FoldersPaths.outputFARP + "", all, nutils.hg, nutils.fpl);
    }

    /**
     *
     * @param all
     * @throws FileNotFoundException
     * @throws IOException
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
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void adenocarcinoma_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlAd, FoldersPaths.fpAd, dbutils.hugo);
        nutils.getLongPathsforListBF(all, dbutils.hugo, 5, FoldersPaths.foundDFAd, "AD");
        putils.rank_pathways(FoldersPaths.foundDFAd, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void adenocarcinoma_part2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.calculate_connectivity_for_pathways(all, FoldersPaths.foundDFAd + "_top_ranked", FoldersPaths.foundDFAd + "_top_ranked_with_connectivity");
        //
        outils.get_centrality_paths_wo_redundancy(all, FoldersPaths.foundDFAd + "_top_ranked", FoldersPaths.foundDFAd + "_top_ranked__unique_overrepr", 2, 4);
        outils.put_names_on_common_parts_wo_redundancy(FoldersPaths.foundDFAd + "_top_ranked__unique_overrepr", FoldersPaths.foundDFAd + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.put_connectivity_on_paths_wo_redundancy(FoldersPaths.foundDFAd + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFAd + "_top_ranked_with_connectivity", FoldersPaths.foundDFAd + "_top_ranked__unique_overrepr_w_names_pthconn");
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param rand
     * @throws IOException
     */
    public void adenocarcinoma_random(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlAd, FoldersPaths.fpAd, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list=new ArrayList(nutils.hg.keySet());
        rand.build_paths_for_random_hit_lists(all, nutils, putils, outils, dbutils, FoldersPaths.Adrandom_path, "AdR", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpAd);
        rand.calculate_p_values(all, nutils, putils, outils, dbutils, FoldersPaths.Adrandom_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFAd + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.Adrandom_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void adenocarcinoma_part1_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlAd, FoldersPaths.fpAd_ppi, dbutils.hugo);
        nutils.getLongPathsforListBF_ppi(all, dbutils.hugo, 5, FoldersPaths.foundDFAd_ppi, "Ad");
        putils.rank_pathways(FoldersPaths.foundDFAd_ppi, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void adenocarcinoma_part2_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.calculate_connectivity_for_pathways(all, FoldersPaths.foundDFAd_ppi + "_top_ranked", FoldersPaths.foundDFAd_ppi + "_top_ranked_with_connectivity");
        outils.get_centrality_paths_wo_redundancy_ppi(all, FoldersPaths.foundDFAd_ppi + "_top_ranked", FoldersPaths.foundDFAd_ppi + "_top_ranked__unique_overrepr", 2, 4);
        outils.put_names_on_common_parts_wo_redundancy(FoldersPaths.foundDFAd_ppi + "_top_ranked__unique_overrepr", FoldersPaths.foundDFAd_ppi + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.put_connectivity_on_paths_wo_redundancy(FoldersPaths.foundDFAd_ppi + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFAd_ppi + "_top_ranked_with_connectivity", FoldersPaths.foundDFAd_ppi + "_top_ranked__unique_overrepr_w_names_pthconn");
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param rand
     * @throws IOException
     */
    public void adenocarcinoma_random_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlAd, FoldersPaths.fpAd_ppi, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list=new ArrayList(nutils.hg.keySet());
        rand.build_paths_for_random_hit_lists(all, nutils, putils, outils, dbutils, FoldersPaths.Adrandom_path_ppi, "AdR", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpAd_ppi);
        rand.calculate_p_values_ppi(all, nutils, putils, outils, dbutils, FoldersPaths.Adrandom_path_ppi + "foundDF.txt_top_ranked", FoldersPaths.foundDFAd_ppi + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.Adrandom_path_ppi + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void mirna63Sys_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        nutils.getLongPathsforListBF(all, dbutils.hugo, 6, FoldersPaths.foundDF63, "MR");
        putils.rank_pathways(FoldersPaths.foundDF63, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void mirna63Sys_part2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.calculate_connectivity_for_pathways(all, FoldersPaths.foundDF63 + "_top_ranked", FoldersPaths.foundDF63 + "_top_ranked_with_connectivity");
        //
        outils.get_centrality_paths_wo_redundancy(all, FoldersPaths.foundDF63 + "_top_ranked", FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr", 1, 4);
        outils.put_names_on_common_parts_wo_redundancy(FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr", FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.put_connectivity_on_paths_wo_redundancy(FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDF63 + "_top_ranked_with_connectivity", FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr_w_names_pthconn");
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param rand
     * @throws IOException
     */
    public void mirna63Sys_random(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list=new ArrayList(nutils.hg.keySet());
        rand.build_paths_for_random_hit_lists(all, nutils, putils, outils, dbutils, FoldersPaths.M63random_path, "63R", 5000, hit_list, dbutils.hugo_by_id, FoldersPaths.fp63);
        rand.calculate_p_values(all, nutils, putils, outils, dbutils, FoldersPaths.M63random_path + "foundDF.txt_top_ranked", FoldersPaths.foundDF63 + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.M63random_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void mirna63Sys_part3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        putils.create_cyto_for_list(FoldersPaths.foundDF63 + "_top_ranked", FoldersPaths.outputF63 + "networks\\paths.txt", FoldersPaths.outputF63 + "networks\\", all, nutils.hg, nutils.fpl);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public static void mirna63Sys_part3_3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        putils.filterPathways_by_length(FoldersPaths.foundDF63, "_max_4", 0, 4);
        outils.find_overreprs(3, 2, FoldersPaths.foundDF63 + "_max_4", FoldersPaths.foundDF63 + "_max_4_overrepr");
        putils.create_cyto_for_paths(FoldersPaths.foundDF63 + "_max_4_overrepr", all, nutils.hg, nutils.fpl);
        putils.filterPathways(FoldersPaths.foundDF63, "_mir_in_the_middle", all, nutils.hg, nutils.fpl, "all", "", "");
        outils.find_overreprs(4, 3, FoldersPaths.foundDF63 + "_mir_in_the_middle", FoldersPaths.foundDF63 + "_mir_in_the_middle_overrepr");
        putils.create_cyto_for_paths(FoldersPaths.foundDF63 + "_mir_in_the_middle_overrepr", all, nutils.hg, nutils.fpl);
        putils.calcmiRNAs(FoldersPaths.foundDF63, FoldersPaths.foundDF63_mirnas, dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "");
        putils.calcmiRNAs(FoldersPaths.foundDF63 + "_mir_in_the_middle", FoldersPaths.foundDF63_mirnas + "_mir_in_the_middle", dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "countonlymiddle");
        String[] mirnas = {"hsa-mir-125b", "hsa-mir-17", "hsa-let-7a"};
        for (String mirna : mirnas) {
            putils.filterPathways_2(FoldersPaths.foundDF63, "_" + mirna.substring(8), all, nutils.hg, nutils.fpl, mirna, 0, 4, 4, 6);
            outils.find_overreprs(5, 3, FoldersPaths.foundDF63 + "_" + mirna.substring(8), FoldersPaths.foundDF63 + "_" + mirna.substring(8) + "_overrepr");
            putils.create_cyto_for_paths(FoldersPaths.foundDF63 + "_" + mirna.substring(8) + "_overrepr", all, nutils.hg, nutils.fpl);
        }
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public static void mirna63Sys_part4(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        outils.filter_paths(FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names", "_min_occ_8", all, nutils.hg, nutils.fpl, 8, "occ");
        outils.filter_paths(FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names_min_occ_8", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.calcmiRNAs_overrep(FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names_min_occ_8_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 8, "");
        outils.filter_paths(FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.calcmiRNAs_overrep(FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 0, "");
        //
        outils.filter_paths(FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names", "_min_occ_20", all, nutils.hg, nutils.fpl, 20, "occ");
        outils.filter_paths(FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_min_occ_20", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.calcmiRNAs_overrep(FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_min_occ_20_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 20, "");
        outils.filter_paths(FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.calcmiRNAs_overrep(FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 0, "");
        putils.create_cyto_for_paths(FoldersPaths.foundDF63 + "_top_ranked_overrepr", all, nutils.hg, nutils.fpl);
    }

 

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void case2_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);
        nutils.getLongPathsforListBF(all, dbutils.hugo, 6, FoldersPaths.foundDFC2, "2C");
        putils.rank_pathways(FoldersPaths.foundDFC2, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void case2_part2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.calculate_connectivity_for_pathways(all, FoldersPaths.foundDFC2 + "_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked_with_connectivity");
        //
        outils.get_centrality_paths_wo_redundancy(all, FoldersPaths.foundDFC2 + "_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr", 1, 4);
        outils.put_names_on_common_parts_wo_redundancy(FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.put_connectivity_on_paths_wo_redundancy(FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFC2 + "_top_ranked_with_connectivity", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names_pthconn");
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void case2_part3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        putils.filterPathways(FoldersPaths.foundDFC2 + "_not_flt", FoldersPaths.foundDFC2, all, nutils.hg, nutils.fpl, "mir_middle", "", "");
        putils.rank_pathways(FoldersPaths.foundDFC2, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
        putils.calculate_connectivity_for_pathways(all, FoldersPaths.foundDFC2 + "_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked_with_connectivity");
        //
        outils.get_centrality_paths_wo_redundancy(all, FoldersPaths.foundDFC2 + "_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr", 1, 4);
        outils.put_names_on_common_parts_wo_redundancy(FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names", all, dbutils.hugo_by_id);
        outils.put_connectivity_on_paths_wo_redundancy(FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names", FoldersPaths.foundDFC2 + "_top_ranked_with_connectivity", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names_pthconn");
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public void case2_part4(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);
        putils.compare_two_overreps(all, dbutils.hugo_by_id, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked__unique_overrepr_w_names_pthconn_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked__unique_overrepr_w_names_pthconn_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\c2_63.txt");
        putils.compare_two_hubs(all, dbutils.hugo_by_id, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked_hubs_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked_hubs_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\c2_63_hubs.txt");
        putils.compare_overreps_with_hitlist(all, nutils, dbutils.hugo_by_id, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked__unique_overrepr_w_names_pthconn_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\c2_hits_63_overreps.txt");
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        putils.compare_overreps_with_hitlist(all, nutils, dbutils.hugo_by_id, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked__unique_overrepr_w_names_pthconn_pvalues", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\63_hits_c2_overreps.txt");
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param rand
     * @throws IOException
     */
    public void case2_random(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);
        int initial_hitlist_size = nutils.hg.size();
        ArrayList<String> hit_list=new ArrayList(nutils.hg.keySet());
        rand.build_paths_for_random_hit_lists(all, nutils, putils, outils, dbutils, FoldersPaths.C2random_path, "C2R", 10000, hit_list, dbutils.hugo_by_id, FoldersPaths.fpC2);
        rand.calculate_p_values(all, nutils, putils, outils, dbutils, FoldersPaths.C2random_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked__unique_overrepr_w_names_pthconn", FoldersPaths.C2random_path + "RandomHitGenes_lists.txt", initial_hitlist_size, 6, dbutils.hugo_by_id);
        rand.calculate_p_values_hubs(all, nutils, putils, outils, dbutils, FoldersPaths.C2random_path + "foundDF.txt_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked_hubs", FoldersPaths.C2random_path + "RandomHitGenes_lists.txt", initial_hitlist_size, dbutils.hugo_by_id);
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param rand
     * @throws IOException
     */
    public void case2_random2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, RandomManager rand) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);

        rand.test_on_random_networks(all, nutils, putils, outils, dbutils, FoldersPaths.C2random_path2 + "foundDF.txt", 1, FoldersPaths.foundDFC2 + "_top_ranked_hubs", "");
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public static void case2_part2_2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);

        outils.filter_paths(FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names", "_min_occ_10", all, nutils.hg, nutils.fpl, 10, "occ");
        outils.filter_paths(FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names_min_occ_10", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.calcmiRNAs_overrep(FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names_min_occ_10_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 10, "");
        outils.filter_paths(FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.calcmiRNAs_overrep(FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 0, "");
        //
        outils.filter_paths(FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr_w_names", "_min_occ_40", all, nutils.hg, nutils.fpl, 40, "occ");
        outils.filter_paths(FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_min_occ_40", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.calcmiRNAs_overrep(FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_min_occ_40_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 40, "");
        outils.filter_paths(FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr_w_names", "_mir_middle", all, nutils.hg, nutils.fpl, 0, "mirna_midlle");
        putils.calcmiRNAs_overrep(FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_mir_middle", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 0, "");
        //
        putils.create_cyto_for_paths(FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr", all, nutils.hg, nutils.fpl);
        putils.create_cyto_for_paths(FoldersPaths.foundDFC2 + "_top_ranked_overrepr", all, nutils.hg, nutils.fpl);
        putils.calcmiRNAs(FoldersPaths.foundDFC2, FoldersPaths.foundDFC2_mirnas, dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "");

    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @throws IOException
     */
    public static void compare_63_c2(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
        /*        nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hlC2, FoldersPaths.fpC2, dbutils.hugo);        
         putils.compare_two_top_ranked(all, dbutils.hugo_by_id,FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names_min_occ_8", FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names_min_occ_10", FoldersPaths.foundDFC2 + "_compare_63_C2_top","63C2");
         putils.compare_two_top_ranked(all, dbutils.hugo_by_id, FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_min_occ_20", FoldersPaths.foundDFC2 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_min_occ_40", FoldersPaths.foundDFC2 + "_compare_63_C2_top_ranked_ppi_min_plus_1","63C2_1");
         putils.find_paths_of_lenght1(FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names", FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names_singles_max");
         putils.find_paths_of_lenght1(FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names", FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names_singles_max");
         putils.create_cyto_for_comparisson(FoldersPaths.foundDFC2 + "_compare_63_C2_top", FoldersPaths.foundDF63 + "_top_ranked", FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names_singles_max", FoldersPaths.foundDFC2 + "_top_ranked", FoldersPaths.foundDFC2 + "_top_ranked_overrepr_w_names_singles_max", FoldersPaths.foundDFC2 + "_compare_63_C2_top_cyto", all);
        
         putils.compare_two_top_ranked(all, dbutils.hugo_by_id, FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names_min_occ_8", FoldersPaths.foundDFC2e + "_top_ranked_overrepr_w_names_min_occ_10", FoldersPaths.foundDFC2e + "_compare_63_C2e_top", "63C2E");
         putils.compare_two_top_ranked(all, dbutils.hugo_by_id, FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names_min_occ_20", FoldersPaths.foundDFC2e + "_top_ranked_ppi_min_plus_1_overrepr_w_names_min_occ_40", FoldersPaths.foundDFC2e + "_compare_63_C2e_top_ranked_ppi_min_plus_1", "63C2e_1");
         putils.find_paths_of_lenght1(FoldersPaths.foundDFC2e + "_top_ranked_overrepr_w_names", FoldersPaths.foundDFC2e + "_top_ranked_overrepr_w_names_singles_max");
         putils.create_cyto_for_comparisson(FoldersPaths.foundDFC2e + "_compare_63_C2e_top", FoldersPaths.foundDF63 + "_top_ranked", FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names_singles_max", FoldersPaths.foundDFC2e + "_top_ranked", FoldersPaths.foundDFC2e + "_top_ranked_overrepr_w_names_singles_max", FoldersPaths.foundDFC2e + "_compare_63_C2e_top_cyto", all);
         */
        /*    putils.compare_two_lists(FoldersPaths.outputFC2E + "genes_top.txt", FoldersPaths.outputF63 + "genes_top.txt", FoldersPaths.outputF63 + "genes_top_C2e_63.txt");
         putils.compare_two_lists(FoldersPaths.outputFC2E + "genes_top_plus_1.txt", FoldersPaths.outputF63 + "genes_top_plus_1.txt", FoldersPaths.outputF63 + "genes_top_plus_1_C2e_63.txt");

         putils.compare_two_lists(FoldersPaths.outputFC2 + "genes_top.txt", FoldersPaths.outputF63 + "genes_top.txt", FoldersPaths.outputF63 + "genes_top_C2_63.txt");
         putils.compare_two_lists(FoldersPaths.outputFC2 + "genes_top_plus_1.txt", FoldersPaths.outputF63 + "genes_top_plus_1.txt", FoldersPaths.outputF63 + "genes_top_plus_1_C2_63.txt");

         putils.compare_two_lists(FoldersPaths.outputFC2E + "reactions_top.txt", FoldersPaths.outputF63 + "reactions_top.txt", FoldersPaths.outputF63 + "reactions_top_C2e_63.txt");
         putils.compare_two_lists(FoldersPaths.outputFC2E + "reactions_top_plus_1.txt", FoldersPaths.outputF63 + "reactions_top_plus_1.txt", FoldersPaths.outputF63 + "reactions_top_plus_1_C2e_63.txt");

         putils.compare_two_lists(FoldersPaths.outputFC2 + "reactions_top.txt", FoldersPaths.outputF63 + "reactions_top.txt", FoldersPaths.outputF63 + "reactions_top_C2_63.txt");
         putils.compare_two_lists(FoldersPaths.outputFC2 + "reactions_top_plus_1.txt", FoldersPaths.outputF63 + "reactions_top_plus_1.txt", FoldersPaths.outputF63 + "reactions_top_plus_1_C2_63.txt");

         putils.compare_two_lists(FoldersPaths.outputFC2 + "genes_top.txt", FoldersPaths.outputFC2E + "genes_top.txt", FoldersPaths.outputFC2E + "genes_top_C2_C2e.txt");
         putils.compare_two_lists(FoldersPaths.outputFC2 + "genes_top_plus_1.txt", FoldersPaths.outputFC2E + "genes_top_plus_1.txt", FoldersPaths.outputFC2E + "genes_top_plus_1_C2_C2e.txt");

         putils.compare_two_lists(FoldersPaths.outputFC2E + "reactions_top.txt", FoldersPaths.outputFC2 + "reactions_top.txt", FoldersPaths.outputFC2E + "reactions_top_C2_C2e.txt");
         putils.compare_two_lists(FoldersPaths.outputFC2E + "reactions_top_plus_1.txt", FoldersPaths.outputFC2 + "reactions_top_plus_1.txt", FoldersPaths.outputFC2E + "reactions_top_plus_1_C2_C2e.txt");
         */
    }

}

/*
 public static void muscleDifSys_part1(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
 nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hl, FoldersPaths.fp, dbutils.hugo);
 //      nutils.getLongPathsforListDF(all, dbutils.hugo, 6, FoldersPaths.foundDF, "HP");
 putils.rank_pathways(FoldersPaths.foundDF, "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);
 putils.rank_pathways(FoldersPaths.foundDF, "_top_ranked_ppi_min_plus_1", all, nutils.hg, nutils.fpl, 1, 1, 0, 0);
 outils.find_overreprs(4, 2, FoldersPaths.foundDF, FoldersPaths.foundDF + "_overrepr");
 outils.find_overreprs(4, 2, FoldersPaths.foundDF + "_top_ranked", FoldersPaths.foundDF + "_top_ranked_overrepr");
 outils.find_overreprs(4, 2, FoldersPaths.foundDF + "_top_ranked_ppi_min_plus_1", FoldersPaths.foundDF + "_top_ranked_ppi_min_plus_1_overrepr");
 outils.put_gene_symbols_for_paths(FoldersPaths.foundDF + "_overrepr", FoldersPaths.foundDF + "_overrepr_w_names", all, dbutils.hugo_by_id);
 outils.put_gene_symbols_for_paths(FoldersPaths.foundDF + "_top_ranked_overrepr", FoldersPaths.foundDF + "_top_ranked_overrepr_w_names", all, dbutils.hugo_by_id);
 outils.put_gene_symbols_for_paths(FoldersPaths.foundDF + "_top_ranked_ppi_min_plus_1_overrepr", FoldersPaths.foundDF + "_top_ranked_ppi_min_plus_1_overrepr_w_names", all, dbutils.hugo_by_id);
 }

 public static void muscleDifSys_part3(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) throws IOException {
 nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hl, FoldersPaths.fp, dbutils.hugo);
 putils.calcmiRNAs(FoldersPaths.foundDF, FoldersPaths.foundDF_mirnas, dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "");
 putils.calcmiRNAs(FoldersPaths.foundDF, FoldersPaths.foundDF_mirnas, dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "countonlymiddle");
 putils.calcmiRNAs(FoldersPaths.foundDF + "_top_ranked", FoldersPaths.foundDF_mirnas + "_top_ranked", dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "");
 putils.calcmiRNAs(FoldersPaths.foundDF + "_top_ranked_ppi_min_plus_1", FoldersPaths.foundDF_mirnas + "_top_ranked_ppi_min_plus_1", dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2, "");
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
 nutils.loadHitGenes_and_FinalPlayers(FoldersPaths.hl, FoldersPaths.fp, dbutils.hugo);
 putils.calcmiRNAs_overrep(FoldersPaths.foundDF + "_overrepr_w_names", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 3, "");
 putils.calcmiRNAs_overrep(FoldersPaths.foundDF + "_top_ranked_overrepr_w_names", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 3, "");
 putils.calcmiRNAs_overrep(FoldersPaths.foundDF + "_top_ranked_ppi_min_plus_1_overrepr_w_names", "", dbutils.mirtarbase, dbutils.transmir, nutils.hg, nutils.fpl, 0, 3, "");
 }

 */
