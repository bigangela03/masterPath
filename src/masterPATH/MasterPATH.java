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
         all0_hippie_h.integrateNetworks(nutils, putils, outils, dbutils);
         all0_hippie_h.saveNetworktofile(FoldersPaths.dbpath + "nodes_hippie_h", FoldersPaths.dbpath + "int_hippie_h");
         all0_hippie_h.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hippie_h", FoldersPaths.dbpath + "int_hippie_h");
         nutils.build_map_of_neighbors(all0_hippie_h);
       
         //PPI network with interactions for which confidence score is more 0
         all0_hippie_m0.integrateNetworks(nutils, putils, outils, dbutils);
         all0_hippie_m0.saveNetworktofile(FoldersPaths.dbpath + "nodes_hippie_m0", FoldersPaths.dbpath + "int_hippie_m0");
         all0_hippie_m0.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hippie_m0", FoldersPaths.dbpath + "int_hippie_m0");
         nutils.build_map_of_neighbors(all0_hippie_m0);
         
         //PPI network with all available interactions
         all0_hippie_a.integrateNetworks(nutils, putils, outils, dbutils);
         all0_hippie_a.saveNetworktofile(FoldersPaths.dbpath + "nodes_hippie_a", FoldersPaths.dbpath + "int_hippie_a");
         all0_hippie_a.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hippie_a", FoldersPaths.dbpath + "int_hippie_a");
         nutils.build_map_of_neighbors(all0_hippie_a);
      
         //PPI network with high + medium confidence interactions        
         all0_hippie_meq.integrateNetworks(nutils, putils, outils, dbutils);
         all0_hippie_meq.saveNetworktofile(FoldersPaths.dbpath + "nodes_hippie_meq", FoldersPaths.dbpath + "int_hippie_meq");
         all0_hippie_meq.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hippie_meq", FoldersPaths.dbpath + "int_hippie_meq");
         nutils.build_map_of_neighbors(all0_hippie_meq);
    
         //PPI network with only HPRD interactions        
         all0_hprd.integrateNetworks(nutils, putils, outils, dbutils);
         all0_hprd.saveNetworktofile(FoldersPaths.dbpath + "nodes_hprd", FoldersPaths.dbpath + "int_hprd");
         all0_hprd.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hprd", FoldersPaths.dbpath + "int_hprd");
         nutils.build_map_of_neighbors(all0_hprd);
     */
         //Integrated network        
         //all2.integrateNetworks(nutils, putils, outils, dbutils);
         //all2.saveNetworktofile(FoldersPaths.dbpath + "all_nodes_all_hippie_high_ogg1", FoldersPaths.dbpath + "all_int_all_hippie_high_ogg1");
         //all2.loadNetworkfromfile(FoldersPaths.dbpath + "all_nodes_all_hippie_high_ogg1", FoldersPaths.dbpath + "all_int_all_hippie_high_ogg1");
         //nutils.build_map_of_neighbors(all2);
       
         //In vivo ppi network        
         //all_invivo.integrateNetworks(nutils, putils, outils, dbutils);
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
        */
         //putils.create_description_for_pathways(FoldersPaths.foundDFOG + "_top_ranked", FoldersPaths.outputFOG + "\\networks\\paths.txt", FoldersPaths.outputFOG + "\\networks\\", all2, nutils.hg,dbutils.hugo_by_id);
         
         /*
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
     

      
         
      */
         //Human muscle differentiation process
        // miRNA loss-of-function screening
        // wrap.mirna63Sys_part1(all2, nutils, putils, outils, dbutils);
        // wrap.mirna63Sys_part2(all2, nutils, putils, outils, dbutils);
        // wrap.mirna63Sys_part3(all2, nutils, putils, outils, dbutils);
        // wrap.mirna63Sys_random(all2, nutils, putils, outils, dbutils,rand);         

/*
        //transcriptomic profiling
         wrap.transcriptomic_part1(all2, nutils, putils, outils, dbutils);
         wrap.transcriptomic_part2(all2, nutils, putils, outils, dbutils);
         wrap.transcriptomic_part3(all2, nutils, putils, outils, dbutils);
         outils.calculate_centrality_scores_for_nodes(all2, dbutils.hugo_by_id, nutils, "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked", "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\foundDF_18_01.txt_top_ranked_hubs");
         wrap.transcriptomic_part4(all2, nutils, putils, outils, dbutils);
         wrap.transcriptomic_random(all2, nutils, putils, outils, dbutils, rand);
         rand.create_random_network(all2, nutils, putils, outils, dbutils);
         wrap.transcriptomic_random2(all2, nutils, putils, outils, dbutils, rand);
     
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

