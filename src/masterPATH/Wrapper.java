package masterPATH;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

/**
 * Wrapper class contains wrappers implementing the whole pipeline
 *
 * @author Natalia Rubanova
 */
public class Wrapper {

    DBManager dbutils;
    NetworkManager nutils = new NetworkManager();
    PathwayManager putils = new PathwayManager();
    CentralityManager outils = new CentralityManager();
    NetworkTopology utils = new NetworkTopology();
    //Wrapper wrap = new Wrapper();
    RandomManager rand = new RandomManager();

    public Wrapper(String hugo_path, String hugo_entrez_path, String entrz_wo_hugo, String hprd_bioDB, String tfacts_crosstable) throws IOException {
        this.dbutils = new DBManager(hugo_path, hugo_entrez_path, entrz_wo_hugo, hprd_bioDB, tfacts_crosstable);
    }

    public Network load_network(
            String file_nodes,
            String file_interaction) throws IOException {
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
            int min_centrality,
            String folder_for_random_paths,
            String prefix_for_random_paths,
            int number_of_permutations,
            int threads) throws IOException, InterruptedException {

        nutils.load_hitlist_and_finalimpl(nw, file_hitlist, file_fimplementers, dbutils.hugo);

        nutils.find_pathway_for_list_BF_algorithm(
                nw,
                dbutils.hugo,
                max_length_for_shortest_path,
                file_output,
                prefix);

        putils.find_the_shortest_paths(
                file_output,
                "_shortest_paths",
                nw,
                nutils.hg,
                nutils.fpl,
                0, 0, 0, 0);

        putils.add_connectivity_to_pathways(
                nw,
                file_output + "_shortest_paths",
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

        ExecutorService es = Executors.newCachedThreadPool();

        NetworkManager[] nutilsrand = new NetworkManager[threads];

        for (int i = 1; i <= threads; i++) {
            int perm_per_thread = 0;
            if (i != threads) {
                perm_per_thread = number_of_permutations / threads;
            } else {
                perm_per_thread = number_of_permutations - (threads - 1) * (number_of_permutations / threads);
            }

            int thread = i;

//            rand.permute_phenotype_label(
//                nw, 
//                nutils, 
//                putils, 
//                outils, 
//                dbutils,
//                folder_for_random_paths, 
//                thread,
//                prefix_for_random_paths + "_" + Integer.toString(perm_per_thread) , 
//                number_of_permutations,
//                hit_list, 
//                dbutils.hugo_by_id, 
//                file_fimplementers);
            // System.out.println("Thread " + thread);
            nutilsrand[i - 1] = new NetworkManager();

            es.execute(new RunnablePermuteLabel(
                    nw,
                    nutilsrand[i - 1],
                    putils,
                    outils,
                    dbutils,
                    folder_for_random_paths,
                    thread,
                    prefix_for_random_paths,
                    number_of_permutations,
                    hit_list,
                    dbutils.hugo_by_id,
                    file_fimplementers,
                    max_length_for_shortest_path
            ));

//            RunnablePermuteLabel myRunnable = new RunnablePermuteLabel(                
//                nw, 
//                nutils, 
//                putils, 
//                outils, 
//                dbutils,
//                folder_for_random_paths, 
//                thread,
//                prefix_for_random_paths + "_" + Integer.toString(perm_per_thread) , 
//                number_of_permutations,
//                hit_list, 
//                dbutils.hugo_by_id, 
//                file_fimplementers,
//                max_length_for_shortest_path
//            );
//            
//            Thread t = new Thread(myRunnable);
//            t.start();
        }

        es.shutdown();
        boolean finished = es.awaitTermination(3, TimeUnit.MINUTES);

        // folder_for_random_paths + "random_paths_" + prefix + "_" + thread_prefix + "_shortest"
        // folder_for_random_paths + "RandomHitGenes" + "_" + thread_prefix
        // folder_for_random_paths + "RandomHitGenes_lists" + "_" + thread_prefix
        String random_paths_file = folder_for_random_paths + "random_paths_" + prefix + "R" + "_shortest";
        String random_hit_genes_file = folder_for_random_paths + "RandomHitGenes";
        String random_hit_genes_lists_file = folder_for_random_paths + "RandomHitGenes_lists";
        BufferedWriter random_path = new BufferedWriter(new FileWriter(random_paths_file));
        BufferedWriter random_hit_genes = new BufferedWriter(new FileWriter(random_hit_genes_file));
        BufferedWriter random_hit_genes_lists = new BufferedWriter(new FileWriter(random_hit_genes_lists_file));
        int thread;
        String thread_prefix;
        String s;
        String[] ss;
        BufferedReader random_path_thread;
        BufferedReader random_hit_genes_thread;
        BufferedReader random_hit_genes_lists_thread;
        ArrayList<String> seenCombinations = new ArrayList();
        String tmpComb;
        for (int i = 1; i <= threads; i++) {
            thread = i;
            thread_prefix = Integer.toString(thread);
            random_path_thread = new BufferedReader(new FileReader(folder_for_random_paths + "random_paths_" + prefix + "R" + "_" + thread_prefix + "_shortest"));
            random_hit_genes_thread = new BufferedReader(new FileReader(folder_for_random_paths + "RandomHitGenes" + "_" + thread_prefix));
            random_hit_genes_lists_thread = new BufferedReader(new FileReader(folder_for_random_paths + "RandomHitGenes_lists" + "_" + thread_prefix));

            while ((s = random_path_thread.readLine()) != null) {
                ss = s.split("\t");
                tmpComb = ss[1] + "\t" + ss[ss.length - 1] + "\t" + Integer.toString(ss.length);
                if (!seenCombinations.contains(tmpComb)) {
                    random_path.write(s + "\n");
                    seenCombinations.add(tmpComb);
                }
            }
            while ((s = random_hit_genes_thread.readLine()) != null) {
                random_hit_genes.write(s + "\n");
            }
            while ((s = random_hit_genes_lists_thread.readLine()) != null) {
                random_hit_genes_lists.write(s + "\n");
            }
            random_path_thread.close();
            random_hit_genes_thread.close();
            random_hit_genes_lists_thread.close();
        }
        random_path.close();
        random_hit_genes.close();
        random_hit_genes_lists.close();

        rand.calculate_p_values_phenotype_label_permutation(
                nw,
                nutils,
                putils,
                outils,
                dbutils,
                folder_for_random_paths + "random_paths_" + prefix + "R_shortest",
                file_output + "_paths_centrality",
                folder_for_random_paths + "RandomHitGenes_lists",
                initial_hitlist_size,
                max_len_for_path + 2,
                dbutils.hugo_by_id);

        rand.calculate_p_values_phenotype_label_permutation_nodes(
                nw,
                nutils,
                putils,
                outils,
                dbutils,
                folder_for_random_paths + "random_paths_" + prefix + "R_shortest",
                file_output + "_nodes_centrality",
                folder_for_random_paths + "RandomHitGenes_lists",
                initial_hitlist_size,
                dbutils.hugo_by_id);
        

        rand.adjust_pvalues(
                nw,
                nutils,
                putils,
                outils,
                dbutils,
                file_output + "_paths_centrality"+ "_pvalues",
                file_output + "_paths_centrality" + "_pvaluesadj",
                min_centrality,
                max_len_for_path + 2,
                "paths"
        );

        rand.adjust_pvalues(
                nw,
                nutils,
                putils,
                outils,
                dbutils,
                file_output + "_nodes_centrality" + "_pvalues",
                file_output + "_nodes_centrality" + "_pvaluesadj",
                min_centrality,
                max_len_for_path + 2,
                "nodes");
        
        System.out.println("Done");
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
            int min_centrality,
            String folder_for_random_paths,
            String prefix_for_random_paths,
            int number_of_permutations,
            int threads) throws IOException, InterruptedException {

        nutils.load_hitlist_and_finalimpl(nw, file_hitlist, file_fimplementers, dbutils.hugo);

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

        ExecutorService es = Executors.newCachedThreadPool();

        NetworkManager[] nutilsrand = new NetworkManager[threads];

        for (int i = 1; i <= threads; i++) {
            int perm_per_thread = 0;
            if (i != threads) {
                perm_per_thread = number_of_permutations / threads;
            } else {
                perm_per_thread = number_of_permutations - (threads - 1) * (number_of_permutations / threads);
            }

            int thread = i;

            nutilsrand[i - 1] = new NetworkManager();

            es.execute(new RunnablePermuteLabel(
                    nw,
                    nutilsrand[i - 1],
                    putils,
                    outils,
                    dbutils,
                    folder_for_random_paths,
                    thread,
                    prefix_for_random_paths,
                    number_of_permutations,
                    hit_list,
                    dbutils.hugo_by_id,
                    file_fimplementers,
                    max_length_for_shortest_path
            ));

        }

        es.shutdown();
        boolean finished = es.awaitTermination(3, TimeUnit.MINUTES);

        // folder_for_random_paths + "random_paths_" + prefix + "_" + thread_prefix + "_shortest"
        // folder_for_random_paths + "RandomHitGenes" + "_" + thread_prefix
        // folder_for_random_paths + "RandomHitGenes_lists" + "_" + thread_prefix
        String random_paths_file = folder_for_random_paths + "random_paths_" + prefix + "R" + "_shortest";
        String random_hit_genes_file = folder_for_random_paths + "RandomHitGenes";
        String random_hit_genes_lists_file = folder_for_random_paths + "RandomHitGenes_lists";
        BufferedWriter random_path = new BufferedWriter(new FileWriter(random_paths_file));
        BufferedWriter random_hit_genes = new BufferedWriter(new FileWriter(random_hit_genes_file));
        BufferedWriter random_hit_genes_lists = new BufferedWriter(new FileWriter(random_hit_genes_lists_file));
        int thread;
        String thread_prefix;
        String s;
        String[] ss;
        BufferedReader random_path_thread;
        BufferedReader random_hit_genes_thread;
        BufferedReader random_hit_genes_lists_thread;
        ArrayList<String> seenCombinations = new ArrayList();
        String tmpComb;
        for (int i = 1; i <= threads; i++) {
            thread = i;
            thread_prefix = Integer.toString(thread);
            random_path_thread = new BufferedReader(new FileReader(folder_for_random_paths + "random_paths_" + prefix + "R" + "_" + thread_prefix + "_shortest"));
            random_hit_genes_thread = new BufferedReader(new FileReader(folder_for_random_paths + "RandomHitGenes" + "_" + thread_prefix));
            random_hit_genes_lists_thread = new BufferedReader(new FileReader(folder_for_random_paths + "RandomHitGenes_lists" + "_" + thread_prefix));

            while ((s = random_path_thread.readLine()) != null) {
                ss = s.split("\t");
                tmpComb = ss[1] + "\t" + ss[ss.length - 1] + "\t" + Integer.toString(ss.length);
                if (!seenCombinations.contains(tmpComb)) {
                    random_path.write(s + "\n");
                    seenCombinations.add(tmpComb);
                }
            }
            while ((s = random_hit_genes_thread.readLine()) != null) {
                random_hit_genes.write(s + "\n");
            }
            while ((s = random_hit_genes_lists_thread.readLine()) != null) {
                random_hit_genes_lists.write(s + "\n");
            }
            random_path_thread.close();
            random_hit_genes_thread.close();
            random_hit_genes_lists_thread.close();
        }
        random_path.close();
        random_hit_genes.close();
        random_hit_genes_lists.close();

        rand.calculate_p_values_phenotype_label_permutation(
                nw,
                nutils,
                putils,
                outils,
                dbutils,
                folder_for_random_paths + "random_paths_" + prefix + "R_shortest",
                file_output + "_paths_centrality",
                folder_for_random_paths + "RandomHitGenes_lists",
                initial_hitlist_size,
                max_len_for_path + 2,
                dbutils.hugo_by_id);

//        int initial_hitlist_size = nutils.hg.size();
//        ArrayList<String> hit_list = new ArrayList(nutils.hg.keySet());
//        rand.permute_phenotype_label(
//                nw, 
//                nutils, 
//                putils, 
//                outils, 
//                dbutils,
//                folder_for_random_paths, 
//                prefix_for_random_paths, 
//                number_of_permutations,
//                hit_list, 
//                dbutils.hugo_by_id, 
//                file_fimplementers
//        );
        rand.calculate_p_values_phenotype_label_permutation_ppi(
                nw,
                nutils,
                putils,
                outils,
                dbutils,
                folder_for_random_paths + "random_paths_" + prefix + "R_shortest",
                file_output + "_paths_centrality",
                folder_for_random_paths + "RandomHitGenes_lists",
                initial_hitlist_size,
                max_len_for_path + 2,
                dbutils.hugo_by_id);

        rand.calculate_p_values_phenotype_label_permutation_nodes(
                nw,
                nutils,
                putils,
                outils,
                dbutils,
                folder_for_random_paths + "random_paths_" + prefix + "R_shortest",
                file_output + "_nodes_centrality",
                folder_for_random_paths + "RandomHitGenes_lists",
                initial_hitlist_size,
                dbutils.hugo_by_id);
        
        rand.adjust_pvalues(
                nw,
                nutils,
                putils,
                outils,
                dbutils,
                file_output + "_paths_centrality"+ "_pvalues",
                file_output + "_paths_centrality" + "_pvaluesadj",
                min_centrality,
                max_len_for_path + 2,
                "paths"
        );

        rand.adjust_pvalues(
                nw,
                nutils,
                putils,
                outils,
                dbutils,
                file_output + "_nodes_centrality" + "_pvalues",
                file_output + "_nodes_centrality" + "_pvaluesadj",
                min_centrality,
                max_len_for_path + 2,
                "nodes");
        
        file = new File(file_output + "_nodes_centrality" + "_pvalues");
        file.delete();
        file = new File(file_output + "_paths_centrality"+ "_pvalues");
        file.delete();
        
        System.out.println("Done");

    }

}
