package masterPATH;


import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Main class
 * 
 * @author Natalia Rubanova
 */
public class MasterPATH {

    /**
     * main method
     * @param args
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {

        DBManager dbutils = new DBManager();
        NetworkManager nutils = new NetworkManager();
        PathwayManager putils = new PathwayManager();
        CentralityManager outils = new CentralityManager();
        NetworkTopology utils = new NetworkTopology();
        Wrapper wrap = new Wrapper();
        RandomManager rand = new RandomManager();
        
       
        Network all0_hippie_h = new Network();
        Network all0_hippie_m0 = new Network();
        Network all0_hprd = new Network();
        Network all0_hippie_m = new Network();
        Network all0_hippie_meq = new Network();
        Network all0_hippie_a = new Network();
        Network all2 = new Network();
        Network all_invivo = new Network();
/*
        
         //PPI network with high confidence interactions
         all0_hippie_h.loadNetwork2(nutils, putils, outils, dbutils);
         all0_hippie_h.saveNetworktofile(FoldersPaths.dbpath + "nodes_hippie_h", FoldersPaths.dbpath + "int_hippie_h");
         all0_hippie_h.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hippie_h", FoldersPaths.dbpath + "int_hippie_h");
         nutils.build_map_of_neighbors(all0_hippie_h);
       
         //PPI network with interactions for which confidence score is more 0
         all0_hippie_m0.loadNetwork2(nutils, putils, outils, dbutils);
         all0_hippie_m0.saveNetworktofile(FoldersPaths.dbpath + "nodes_hippie_m0", FoldersPaths.dbpath + "int_hippie_m0");
         all0_hippie_m0.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hippie_m0", FoldersPaths.dbpath + "int_hippie_m0");
         nutils.build_map_of_neighbors(all0_hippie_m0);
         
         //PPI network with all available interactions
         all0_hippie_a.loadNetwork2(nutils, putils, outils, dbutils);
         all0_hippie_a.saveNetworktofile(FoldersPaths.dbpath + "nodes_hippie_a", FoldersPaths.dbpath + "int_hippie_a");
         all0_hippie_a.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hippie_a", FoldersPaths.dbpath + "int_hippie_a");
         nutils.build_map_of_neighbors(all0_hippie_a);
      
         //PPI network with high + medium confidence interactions        
         all0_hippie_meq.loadNetwork2(nutils, putils, outils, dbutils);
         all0_hippie_meq.saveNetworktofile(FoldersPaths.dbpath + "nodes_hippie_meq", FoldersPaths.dbpath + "int_hippie_meq");
         all0_hippie_meq.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hippie_meq", FoldersPaths.dbpath + "int_hippie_meq");
         nutils.build_map_of_neighbors(all0_hippie_meq);
    
         //PPI network with only HPRD interactions        
         all0_hprd.loadNetwork2(nutils, putils, outils, dbutils);
         all0_hprd.saveNetworktofile(FoldersPaths.dbpath + "nodes_hprd", FoldersPaths.dbpath + "int_hprd");
         all0_hprd.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hprd", FoldersPaths.dbpath + "int_hprd");
         nutils.build_map_of_neighbors(all0_hprd);
     */
         //Integrated network        
         //all2.loadNetwork2(nutils, putils, outils, dbutils);
         //all2.saveNetworktofile(FoldersPaths.dbpath + "all_nodes_all_hippie_high_ogg1", FoldersPaths.dbpath + "all_int_all_hippie_high_ogg1");
         //all2.loadNetworkfromfile(FoldersPaths.dbpath + "all_nodes_all_hippie_high_ogg1", FoldersPaths.dbpath + "all_int_all_hippie_high_ogg1");
         //nutils.build_map_of_neighbors(all2);
       
         //In vivo ppi network        
         //all_invivo.loadNetwork2(nutils, putils, outils, dbutils);
         //all_invivo.saveNetworktofile(FoldersPaths.dbpath + "all_nodes_invivo", FoldersPaths.dbpath + "all_int_invivo");
         //all_invivo.loadNetworkfromfile(FoldersPaths.dbpath + "all_nodes_invivo", FoldersPaths.dbpath + "all_int_invivo");
         //nutils.build_map_of_neighbors(all_invivo);
         
         /*
         //DNA repair screening hit genes - OGG1
         wrap.ogg1_part1(all2, nutils, putils, outils, dbutils);
         wrap.ogg1_part2(all2, nutils, putils, outils, dbutils);
         wrap.ogg1_random(all2, nutils, putils, outils, dbutils, rand);
         wrap.ogg1_part3(all2, nutils, putils, outils, dbutils);
         wrap.ogg1_part4(all2, nutils, putils, outils, dbutils);

         //DNA repair screening hit genes - hit genes
         wrap.ogg1_part1_ppi(all0_hippie_h, nutils, putils, outils, dbutils);
         wrap.ogg1_part2_ppi(all0_hippie_h, nutils, putils, outils, dbutils);
         wrap.ogg1_random_ppi(all0_hippie_h, nutils, putils, outils, dbutils, rand);
         wrap.ogg1_part4_ppi(all0_hippie_h, nutils, putils, outils, dbutils);
         
         //Further analysis
         outils.calculate_centrality_scores_for_nodes(all0_hippie_h, dbutils.hugo_by_id, nutils, FoldersPaths.foundDFOG_ppi + "_top_ranked", FoldersPaths.foundDFOG_ppi + "_top_ranked_hubs");
         nutils.load_hitlist_and_finalimpl(FoldersPaths.hlOG, FoldersPaths.fpOG_ppi, dbutils.hugo);
         int initial_hitlist_size = nutils.hg.size();
         rand.calculate_p_values_phenotype_label_permutation_nodes(all0_hippie_h, nutils, putils, outils, dbutils, FoldersPaths.OGrandom_path_ppi + "foundDF.txt_top_ranked", FoldersPaths.foundDFOG_ppi + "_top_ranked_hubs", FoldersPaths.OGrandom_path_ppi + "RandomHitGenes_lists.txt", initial_hitlist_size, dbutils.hugo_by_id);
         wrap.ogg1_random2_ppi(all0_hippie_h, nutils, putils, outils, dbutils, rand);
         nutils.load_hitlist_and_finalimpl(FoldersPaths.hlOG, FoldersPaths.fpOG, dbutils.hugo);
         putils.find_the_shortest_paths(FoldersPaths.foundDFOG, "_top_ranked", all2, nutils.hg, nutils.fpl, 0, 0, 0, 0);
         outils.rank_hit_genes(all2, FoldersPaths.foundDFOG_ppi + "_top_ranked", FoldersPaths.foundDFOG_ppi + "_ranked_hits", dbutils);
     
         //ARP2/3 loss of function screening
         //PPI
         wrap.arp_part1_ppi(all0_hippie_h, nutils, putils, outils, dbutils);
         wrap.arp_part2_ppi(all0_hippie_h, nutils, putils, outils, dbutils);
         wrap.arp_random_ppi(all0_hippie_h, nutils, putils, outils, dbutils, rand); 
         //Integrated network
         wrap.arp_part1_LIMCH1(all2, nutils, putils, outils, dbutils);
         wrap.arp_part1(all2, nutils, putils, outils, dbutils);
         wrap.arp_part3(all2, nutils, putils, outils, dbutils);
         wrap.arp_part2(all2, nutils, putils, outils, dbutils);
         wrap.arp_random(all2, nutils, putils, outils, dbutils, rand);
         wrap.arp_filter(all2); 
      
         //Human adenocarcinoma 
         //Integrated network
         wrap.adenocarcinoma_part1(all2, nutils, putils, outils, dbutils);
         wrap.adenocarcinoma_part2(all2, nutils, putils, outils, dbutils);
         wrap.adenocarcinoma_random(all2, nutils, putils, outils, dbutils, rand);
         //PPI
         wrap.adenocarcinoma_part1_ppi(all0_hippie_h, nutils, putils, outils, dbutils);
         wrap.adenocarcinoma_part2_ppi(all0_hippie_h, nutils, putils, outils, dbutils);
         wrap.adenocarcinoma_random_ppi(all0_hippie_h, nutils, putils, outils, dbutils, rand);
      */
         //Human muscle differentiation process
        // miRNA loss-of-function screening
        // wrap.mirna63Sys_part1(all2, nutils, putils, outils, dbutils);
        // wrap.mirna63Sys_part2(all2, nutils, putils, outils, dbutils);
        // wrap.mirna63Sys_part3(all2, nutils, putils, outils, dbutils);
        // wrap.mirna63Sys_random(all2, nutils, putils, outils, dbutils,rand);         
        // outils.calculate_centrality_scores_for_nodes(all2, dbutils.hugo_by_id, nutils, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked_hubs");
        // nutils.load_hitlist_and_finalimpl(FoldersPaths.hl63, FoldersPaths.fp63, dbutils.hugo);
        // int initial_hitlist_size = nutils.hg.size();
        // rand.calculate_p_values_phenotype_label_permutation_nodes(all2, nutils, putils, outils, dbutils, FoldersPaths.M63random_path + "foundDF.txt_top_ranked", FoldersPaths.foundDF63 + "_top_ranked_hubs", FoldersPaths.M63random_path + "RandomHitGenes_lists.txt", initial_hitlist_size, dbutils.hugo_by_id);
         
        
/*
        //transcriptomic profiling
         wrap.case2_part1(all2, nutils, putils, outils, dbutils);
         wrap.case2_part2(all2, nutils, putils, outils, dbutils);
         wrap.case2_part3(all2, nutils, putils, outils, dbutils);
         outils.calculate_centrality_scores_for_nodes(all2, dbutils.hugo_by_id, nutils, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked_hubs");
         wrap.case2_part4(all2, nutils, putils, outils, dbutils);
         wrap.case2_random(all2, nutils, putils, outils, dbutils, rand);
         rand.create_random_network(all2, nutils, putils, outils, dbutils);
         wrap.case2_random2(all2, nutils, putils, outils, dbutils, rand);
     
         //Networks topological properties
         System.out.println("hprd");
         utils.calculate_average_clustering_coefficient_ppi(all0_hprd);
         utils.calculate_number_of_connected_components_ppi(all0_hprd);
         System.out.println("hippi high");
         utils.calculate_average_clustering_coefficient_ppi(all0_hippie_h);
         utils.calculate_number_of_connected_components_ppi(all0_hippie_h);
         System.out.println("hippi all");
         utils.calculate_average_clustering_coefficient_ppi(all0_hippie_a);
         utils.calculate_number_of_connected_components_ppi(all0_hippie_a);
         System.out.println("hippi meq");
         utils.calculate_average_clustering_coefficient_ppi(all0_hippie_meq);
         utils.calculate_number_of_connected_components_ppi(all0_hippie_meq);
         System.out.println("hprd");
         utils.calculate_diameter_ppi(all0_hprd);
         System.out.println("hippi high");
         utils.calculate_diameter_ppi(all0_hippie_h);
         System.out.println("hippi meq");
         utils.calculate_diameter_ppi(all0_hippie_meq);
         System.out.println("hippi all");
         utils.calculate_diameter_ppi(all0_hippie_a);
         System.out.println("hprd");
         utils.get_degree_distribution_ppi(all0_hprd,"F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\hprd.txt");
         System.out.println("hippi high");
         utils.get_degree_distribution_ppi(all0_hippie_h,"F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\hippie_h.txt");
         System.out.println("hippi meq");
         utils.get_degree_distribution_ppi(all0_hippie_meq,"F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\hippie_meq.txt");        
         System.out.println("hippi all");
         utils.get_degree_distribution_ppi(all0_hippie_a,"F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\hippie_a.txt");        
         System.out.println("all2");
         utils.calculate_number_of_connected_components_direct(all2);
         utils.get_degree_distribution_direct(all2,"F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\all2.txt");
         utils.calculate_average_clustering_coefficient_direct(all2);
         utils.calculate_diameter_direct(all2);
         utils.compare_two_lists();
         utils.get_subnetwork_statistics(all2, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked");
         utils.get_subnetwork_statistics(all2, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\63mirnas\\output\\foundDF_18_01.txt_top_ranked"); 
         utils.get_subnetwork_statistics(all0_hippie_h, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\OGG1\\output\\ppi\\foundDF_18_01.txt_top_ranked");               
         */
    }

}


/*
 //putils.filterPathways(FoldersPaths.foundDFC2, "_mir_all", all, nutils.hg, nutils.fpl, "mir_all", "", "");
 //nutils.find_overreprs(4, 3, FoldersPaths.foundDFC2 + "_mir_all", FoldersPaths.foundDFC2 + "_mir_all_overrepr");
 //putils.create_cyto_for_paths(FoldersPaths.foundDFC2 + "_mir_all_overrepr", all, nutils.hg, nutils.fpl);
 // putils.filterPathways(FoldersPaths.foundDFC2 + "_top_ranked", "_mir_all", all, nutils.hg, nutils.fpl, "mir_all", "", "");

 // nutils.find_overreprs(5, 1, FoldersPaths.foundDF63, FoldersPaths.foundDF63 + "_5_1_overrepr");
 // nutils.find_overreprs(5, 1, FoldersPaths.foundDF63 + "_top_ranked", FoldersPaths.foundDF63 + "_top_ranked_5_1_overrepr");
 // nutils.find_overreprs(5, 1, FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1", FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_5_1_overrepr");


 //foundDF_18_01.txt_mirna_ends_wo_MIRT
 //  putils.compare_two_top_ranked(all, dbutils.hugo_by_id, FoldersPaths.foundDF63 + "_5_1_overrepr", FoldersPaths.foundDFC2 + "_mirna_ends_wo_MIRT", FoldersPaths.foundDFC2 + "_compare_C2_mir_all_w_63_5_1_mirna_ends_wo_MIRT");
 //  putils.compare_two_top_ranked(all, dbutils.hugo_by_id, FoldersPaths.foundDF63 + "_5_1_overrepr", FoldersPaths.foundDFC2 + "_mirna_ends", FoldersPaths.foundDFC2 + "_compare_C2_mir_all_w_63_5_1_mirna_ends");
 //   putils.compare_two_top_ranked(all, dbutils.hugo_by_id, FoldersPaths.foundDF63 + "_top_ranked_5_1_overrepr", FoldersPaths.foundDFC2 + "_mirna_ends_wo_MIRT", FoldersPaths.foundDFC2 + "_compare_C2_mir_all_w_63_5_1_mirna_ends_wo_MIRT");
 //   putils.compare_two_top_ranked(all, dbutils.hugo_by_id, FoldersPaths.foundDF63 + "_top_ranked_5_1_overrepr", FoldersPaths.foundDFC2 + "_mirna_ends", FoldersPaths.foundDFC2 + "_compare_C2_mir_all_w_63_5_1_mirna_ends");
 //   putils.compare_two_top_ranked(all, dbutils.hugo_by_id, FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_5_1_overrepr", FoldersPaths.foundDFC2 + "_mirna_ends_wo_MIRT", FoldersPaths.foundDFC2 + "_compare_C2_mir_all_w_63_ppi_min_plus_1_5_1_mirna_ends_wo_MIRT");
 //   putils.compare_two_top_ranked(all, dbutils.hugo_by_id, FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_5_1_overrepr", FoldersPaths.foundDFC2 + "_mirna_ends", FoldersPaths.foundDFC2 + "_compare_C2_mir_all_w_63_ppi_min_plus_1_5_1_mirna_ends");
 // nutils.put_gene_symbols_for_paths(FoldersPaths.foundDF63 + "_overrepr", FoldersPaths.foundDF63 + "_overrepr_w_names", all, dbutils.hugo_by_id);
 // nutils.put_gene_symbols_for_paths(FoldersPaths.foundDF63 + "_top_ranked_overrepr", FoldersPaths.foundDF63 + "_top_ranked_overrepr_w_names", all, dbutils.hugo_by_id);
 // nutils.put_gene_symbols_for_paths(FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr", FoldersPaths.foundDF63 + "_top_ranked_ppi_min_plus_1_overrepr_w_names", all, dbutils.hugo_by_id);




 // putils.filterPathways_cyto_by_length(FoldersPaths.outputF + "foundDF_28_11.txt_206", "_min_6", 6, 6);

 putils.find_miRNAs_on_pathways(FoldersPaths.foundDF, FoldersPaths.foundDF_mirnas, dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2);
 putils.filterPathways(FoldersPaths.foundDF, "_filter_JAG1mirna", all, nutils.hg, nutils.fpl, "HGNC:6188", "miRNA");
 putils.filterPathways(FoldersPaths.foundDF, "_filter_allmirna", all, nutils.hg, nutils.fpl, "all", "miRNA");
 nutils.find_overreprs(5, 3, FoldersPaths.foundDF + "_filter_JAG1mirna", FoldersPaths.foundDF + "_filter_JAG1mirna" + "_overrepr");
 nutils.find_overreprs(5, 3, FoldersPaths.foundDF + "_filter_allmirna", FoldersPaths.foundDF + "_filter_allmirna" + "_overrepr");
 putils.create_network_for_cytoscape(FoldersPaths.foundDF_overrepr + "_sorted", all, nutils.hg, nutils.fpl);
 putils.create_network_for_cytoscape(FoldersPaths.foundDF + "_filter_JAG1mirna" + "_an", all, nutils.hg, nutils.fpl);
 putils.create_network_for_cytoscape(FoldersPaths.foundDF + "_filter_allmirna" + "_an", all, nutils.hg, nutils.fpl); */

/*
 putils.filterPathways_cyto_by_length(FoldersPaths.outputF + "foundDF_28_11.txt_206", "_max_5", 0, 5);
 nutils.find_overreprs(4, 3, FoldersPaths.outputF + "foundDF_28_11.txt_206_max_5", FoldersPaths.outputF + "foundDF_28_11.txt_206_max_5_overrepr");
 putils.create_cyto_for_paths(FoldersPaths.outputF + "foundDF_28_11.txt_206_max_5_overrepr", all, nutils.hg, nutils.fpl);
        
 putils.filterPathways(FoldersPaths.outputF + "foundDF_28_11.txt_206", "_filter_mirna", all, nutils.hg, nutils.fpl, "all", "", "");
 nutils.find_overreprs(4, 3, FoldersPaths.outputF + "foundDF_28_11.txt_206_filter_mirna", FoldersPaths.outputF + "foundDF_28_11.txt_206_filter_mirna_overrepr");
 putils.create_cyto_for_paths(FoldersPaths.outputF + "foundDF_28_11.txt_206_filter_mirna_overrepr", all, nutils.hg, nutils.fpl);
 */
/*
 String[] mirnas = {"hsa-mir-34a", "hsa-mir-200b", "hsa-mir-27a", "hsa-mir-449a", "hsa-mir-30a", "hsa-mir-148a", "hsa-mir-145", "hsa-mir-206", "hsa-mir-125b", "hsa-mir-125b", "hsa-let-7a"};
 for (String mirna : mirnas) {
 putils.filterPathways(FoldersPaths.foundDF, "_" + mirna.substring(8), all, nutils.hg, nutils.fpl, "mirna", mirna, "miRNA");
 nutils.find_overreprs(5, 3, FoldersPaths.foundDF + "_" + mirna.substring(8), FoldersPaths.foundDF + "_" + mirna.substring(8) + "_overrepr");
 putils.create_cyto_for_paths(FoldersPaths.foundDF + "_" + mirna.substring(8) + "_overrepr", all, nutils.hg, nutils.fpl);
 } */
/*
 putils.filterPathways_cyto_by_length(FoldersPaths.outputF + "foundDF_28_11.txt_206", "_max_5", 0, 5);
 nutils.find_overreprs(4, 3, FoldersPaths.outputF + "foundDF_28_11.txt_206_max_5", FoldersPaths.outputF + "foundDF_28_11.txt_206_max_5_overrepr");
 putils.create_cyto_for_paths(FoldersPaths.outputF + "foundDF_28_11.txt_206_max_5_overrepr", all, nutils.hg, nutils.fpl);
        
 putils.filterPathways(FoldersPaths.outputF + "foundDF_28_11.txt_206", "_filter_mirna", all, nutils.hg, nutils.fpl, "all", "", "");
 nutils.find_overreprs(4, 3, FoldersPaths.outputF + "foundDF_28_11.txt_206_filter_mirna", FoldersPaths.outputF + "foundDF_28_11.txt_206_filter_mirna_overrepr");
 putils.create_cyto_for_paths(FoldersPaths.outputF + "foundDF_28_11.txt_206_filter_mirna_overrepr", all, nutils.hg, nutils.fpl);
 */
/*
 putils.find_miRNAs_on_pathways(FoldersPaths.foundDF, FoldersPaths.foundDF_mirnas, dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2);
 putils.filterPathways(FoldersPaths.foundDF, "_filter_JAG1mirna", all, nutils.hg, nutils.fpl, "HGNC:6188", "miRNA");
 putils.filterPathways(FoldersPaths.foundDF, "_filter_allmirna", all, nutils.hg, nutils.fpl, "all", "miRNA");
 nutils.find_overreprs(5, 3, FoldersPaths.foundDF + "_filter_JAG1mirna", FoldersPaths.foundDF + "_filter_JAG1mirna" + "_overrepr");
 nutils.find_overreprs(5, 3, FoldersPaths.foundDF + "_filter_allmirna", FoldersPaths.foundDF + "_filter_allmirna" + "_overrepr");
 putils.create_network_for_cytoscape(FoldersPaths.foundDF_overrepr + "_sorted", all, nutils.hg, nutils.fpl);
 putils.create_network_for_cytoscape(FoldersPaths.foundDF + "_filter_JAG1mirna" + "_an", all, nutils.hg, nutils.fpl);
 putils.create_network_for_cytoscape(FoldersPaths.foundDF + "_filter_allmirna" + "_an", all, nutils.hg, nutils.fpl); */

/*
 putils.find_miRNAs_on_pathways(FoldersPaths.foundDF, FoldersPaths.foundDF_mirnas, dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2);
 putils.filterPathways(FoldersPaths.foundDF, "_filter_JAG1mirna", all, nutils.hg, nutils.fpl, "HGNC:6188", "miRNA");
 putils.filterPathways(FoldersPaths.foundDF, "_filter_allmirna", all, nutils.hg, nutils.fpl, "all", "miRNA");
 nutils.find_overreprs(5, 3, FoldersPaths.foundDF + "_filter_JAG1mirna", FoldersPaths.foundDF + "_filter_JAG1mirna" + "_overrepr");
 nutils.find_overreprs(5, 3, FoldersPaths.foundDF + "_filter_allmirna", FoldersPaths.foundDF + "_filter_allmirna" + "_overrepr");
 putils.create_network_for_cytoscape(FoldersPaths.foundDF_overrepr + "_sorted", all, nutils.hg, nutils.fpl);
 putils.create_network_for_cytoscape(FoldersPaths.foundDF + "_filter_JAG1mirna" + "_an", all, nutils.hg, nutils.fpl);
 putils.create_network_for_cytoscape(FoldersPaths.foundDF + "_filter_allmirna" + "_an", all, nutils.hg, nutils.fpl); */

/*
 String[] mirnas = {"hsa-mir-34a", "hsa-mir-200b", "hsa-mir-27a", "hsa-mir-449a", "hsa-mir-30a", "hsa-mir-148a", "hsa-mir-145", "hsa-mir-206", "hsa-mir-125b", "hsa-mir-125b", "hsa-let-7a"};
 for (String mirna : mirnas) {
 putils.filterPathways(FoldersPaths.foundDF, "_" + mirna.substring(8), all, nutils.hg, nutils.fpl, "mirna", mirna, "miRNA");
 nutils.find_overreprs(5, 3, FoldersPaths.foundDF + "_" + mirna.substring(8), FoldersPaths.foundDF + "_" + mirna.substring(8) + "_overrepr");
 putils.create_cyto_for_paths(FoldersPaths.foundDF + "_" + mirna.substring(8) + "_overrepr", all, nutils.hg, nutils.fpl);
 } */
/*
 putils.filterPathways_cyto_by_length(FoldersPaths.outputF + "foundDF_28_11.txt_206", "_max_5", 0, 5);
 nutils.find_overreprs(4, 3, FoldersPaths.outputF + "foundDF_28_11.txt_206_max_5", FoldersPaths.outputF + "foundDF_28_11.txt_206_max_5_overrepr");
 putils.create_cyto_for_paths(FoldersPaths.outputF + "foundDF_28_11.txt_206_max_5_overrepr", all, nutils.hg, nutils.fpl);
        
 putils.filterPathways(FoldersPaths.outputF + "foundDF_28_11.txt_206", "_filter_mirna", all, nutils.hg, nutils.fpl, "all", "", "");
 nutils.find_overreprs(4, 3, FoldersPaths.outputF + "foundDF_28_11.txt_206_filter_mirna", FoldersPaths.outputF + "foundDF_28_11.txt_206_filter_mirna_overrepr");
 putils.create_cyto_for_paths(FoldersPaths.outputF + "foundDF_28_11.txt_206_filter_mirna_overrepr", all, nutils.hg, nutils.fpl);
 */
/*
 putils.find_miRNAs_on_pathways(FoldersPaths.foundDF, FoldersPaths.foundDF_mirnas, dbutils.mirtarbase, nutils.hg, nutils.fpl, 6, 2);
 putils.filterPathways(FoldersPaths.foundDF, "_filter_JAG1mirna", all, nutils.hg, nutils.fpl, "HGNC:6188", "miRNA");
 putils.filterPathways(FoldersPaths.foundDF, "_filter_allmirna", all, nutils.hg, nutils.fpl, "all", "miRNA");
 nutils.find_overreprs(5, 3, FoldersPaths.foundDF + "_filter_JAG1mirna", FoldersPaths.foundDF + "_filter_JAG1mirna" + "_overrepr");
 nutils.find_overreprs(5, 3, FoldersPaths.foundDF + "_filter_allmirna", FoldersPaths.foundDF + "_filter_allmirna" + "_overrepr");
 putils.create_network_for_cytoscape(FoldersPaths.foundDF_overrepr + "_sorted", all, nutils.hg, nutils.fpl);
 putils.create_network_for_cytoscape(FoldersPaths.foundDF + "_filter_JAG1mirna" + "_an", all, nutils.hg, nutils.fpl);
 putils.create_network_for_cytoscape(FoldersPaths.foundDF + "_filter_allmirna" + "_an", all, nutils.hg, nutils.fpl); */
/*


 /*String path = "F:\\Dropbox\\_Работа\\Programms\\DATA\\data_for_Merge\\";
 String musclef_com = "F:\\Dropbox\\_Работа\\Programms\\DATA\\Muscle Differentiation\\";
 String musclef = "F:\\Dropbox\\_Работа\\Programms\\DATA\\Muscle Differentiation\\Sep15\\";
 String hl = musclef_com + "HitGenes.txt";
 String fp = musclef_com + "FinalPlayers.txt";
 // String muscle_paths = musclef + "MPaths_22_09.txt";
 String foundDF = musclef + "foundDF_22_09.txt";
 String foundDFan = musclef + "foundDFan_22_09.txt";
 String foundDFmirnas = musclef + "foundDFmirnas_22_09.txt";
 String muscle_nw2 = musclef + "nw2_edt.txt";
 String muscle_nw2_w_connectivities = musclef + "nw2_edt_con.txt"; */

/*
 if (all.nodes.containsKey("pHGNC:12770")) {
 System.out.println("ppprotein");
 System.out.println(all.nodes.get("pHGNC:12770").id_type + "\t " + all.nodes.get("pHGNC:12770").type + "\t " + all.nodes.get("pHGNC:12770").db_flag);
 System.out.println("down");
 for (String s : all.nodes.get("pHGNC:12770").downnbrs.keySet()) {
 System.out.println(s);
 }
 System.out.println("up");

 for (String s : all.nodes.get("pHGNC:12770").upnbrs.keySet()) {
 System.out.println(s);
 }
 System.out.println("rev");

 for (String s : all.nodes.get("pHGNC:12770").revnbrs.keySet()) {
 System.out.println(s);
 }
 }
 if (all.nodes.containsKey("HGNC:12770")) {
 System.out.println("gggene");
 System.out.println(all.nodes.get("HGNC:12770").id_type + "\t " + all.nodes.get("HGNC:12770").type + "\t " + all.nodes.get("HGNC:12770").db_flag);
 System.out.println("down");
 for (String s : all.nodes.get("HGNC:12770").downnbrs.keySet()) {
 System.out.println(s);
 }
 System.out.println("up");

 for (String s : all.nodes.get("HGNC:12770").upnbrs.keySet()) {
 System.out.println(s);
 }
 System.out.println("rev");

 for (String s : all.nodes.get("HGNC:12770").revnbrs.keySet()) {
 System.out.println(s);
 }
 }

 if (all.nodes.containsKey("pHGNC:12770")) {
 System.out.println("ppprotein");
 System.out.println(all.nodes.get("pHGNC:12770").id_type + "\t " + all.nodes.get("pHGNC:12770").type + "\t " + all.nodes.get("pHGNC:12770").db_flag);
 System.out.println("down");
 for (String s : all.nodes.get("pHGNC:12770").downnbrs.keySet()) {
 System.out.println(s);
 }
 System.out.println("up");

 for (String s : all.nodes.get("pHGNC:12770").upnbrs.keySet()) {
 System.out.println(s);
 }
 System.out.println("rev");

 for (String s : all.nodes.get("pHGNC:12770").revnbrs.keySet()) {
 System.out.println(s);
 }
 }
 if (all.nodes.containsKey("HGNC:12770")) {
 System.out.println("gggene");
 System.out.println(all.nodes.get("HGNC:12770").id_type + "\t " + all.nodes.get("HGNC:12770").type + "\t " + all.nodes.get("HGNC:12770").db_flag);
 System.out.println("down");
 for (String s : all.nodes.get("HGNC:12770").downnbrs.keySet()) {
 System.out.println(s);
 }
 System.out.println("up");

 for (String s : all.nodes.get("HGNC:12770").upnbrs.keySet()) {
 System.out.println(s);
 }
 System.out.println("rev");

 for (String s : all.nodes.get("HGNC:12770").revnbrs.keySet()) {
 System.out.println(s);
 }
 }



 NetworkManager nutils = new NetworkManager();
 nutils.load_hitlist_and_finalimpl(hl, fp, dbutils.hugo);
 List<PathUnit> npaths = new ArrayList();
 List<PathUnit> npathsBF = new ArrayList();
 List<PathUnit> npathsDF = new ArrayList();
 // nutils.getPathsforList(all, hl, fp, dbutils.hugo, 3, muscle_paths);
 // nutils.find_pathway_for_list_BF_algorithm(all, hl, fp, dbutils.hugo, 5, musclef);
 // nutils.getConnectivity(all, muscle_nw2, muscle_nw2_w_connectivities);
 // npaths = nutils.getPath(all, "pHGNC:7611", "pHGNC:6996", 3);
 // npathsBF=nutils.find_pathway_BF_algorithm(all, "HGNC:7611", "pHGNC:6996", 3);
 // npathsDF = nutils.getLongPathDF(all, "pHGNC:7611", "pHGNC:6996", 5);
 //!!!    nutils.getLongPathsforListDF(all, hl, fp, dbutils.hugo, 5, foundDF);
 System.out.println("npaths size " + npaths.size() + " npathsDF size " + npathsDF.size());               
 System.out.println("npaths");
 for (PathUnit np : npaths) {
 for (int j = 0; j < np.path.size(); j++) {
 System.out.print(np.path.get(j) + " ");
 }
 System.out.print("\n");
 }
 System.out.println("npathsDF");
 for (PathUnit np : npathsDF) {
 for (int j = 0; j < np.path.size(); j++) {
 System.out.print(np.path.get(j) + " ");
 }
 System.out.print("\n");
 }
         
 //System.out.println(npaths.size());


 for (int i = 0; i < npaths.size(); i++) {
 for (int j = 0; j < npaths.get(i).size(); j++) {
 System.out.print(all.interactions.get(npaths.get(i).get(j)).sourcedbentry.get(0)[1] + "\t" + all.interactions.get(npaths.get(i).get(j)).sourcedbentry.get(0)[1] + "-" +all.interactions.get(npaths.get(i).get(j)).sourcedbentry.get(0)[4] +"\t-------");
 }
 System.out.print("\n");
 } 

 */
