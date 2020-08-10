package masterpath.masterpath;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

/**
 * RandomManager contains methods to calculate p-values
 *
 * @author Natalia Rubanova
 */
public class RandomManager {

    /**
     * Permute phenotype label
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     * @param folder Folder to store information about pathways
     * @param prefix Prefix for pathways ids
     * @param count Number of shuffling steps
     * @param hit_list Hit List
     * @param hugo_by_id HGNC ids map
     * @param fplayers Final implementers
     */
    public void permute_phenotype_label(
            Network all,
            NetworkManager nutils,
            PathwayManager putils,
            CentralityManager outils,
            DBManager dbutils,
            String folder,
            int thread,
            String prefix,
            int count,
            List<String> hit_list,
            Map<String, String[]> hugo_by_id,
            String fplayers,
            int maxLength
    ) throws IOException {

        System.out.println("+++++++++++++build_paths_for_random_hit_lists++++++++++++");
        java.util.Random randomGenerator = new java.util.Random();
        List<String> all_hugo_ids = new ArrayList();
        List<String[]> all_mirnas_ids = new ArrayList();
        List<String[]> all_chemicals_ids = new ArrayList();
        int count_list_size = hit_list.size();
        int randomInt, degree;
        int max_degree = 0;
        int bin = 0;
        int index = 0;
        int number_of_proteins = 0;
        String thread_prefix = Integer.toString(thread);
        for (Node n : all.nodes.values()) {
            if (n.type.equals("protein")) {
                degree = n.downnbrs.size() + n.upnbrs.size() + n.revnbrs.size();
                if (degree > 0) {
                    number_of_proteins++;
                }
            }
        }
        int[] degree_distribution = new int[number_of_proteins];
        for (Node n : all.nodes.values()) {
            if (n.type.equals("protein")) {
                degree = n.downnbrs.size() + n.upnbrs.size() + n.revnbrs.size();
                if (degree > 0) {
                    degree_distribution[index] = degree;
                    index++;
                    if (max_degree < degree) {
                        max_degree = degree;
                    }
                }
            }
        }

        Arrays.sort(degree_distribution);
        int prev = degree_distribution[0];
        int number = 1;
        for (int i = 1; i < degree_distribution.length; i++) {
            if (prev != degree_distribution[i]) {
                //System.out.print(prev + "-" + number + " ");
                prev = degree_distribution[i];
                number = 1;
            } else {
                number++;
            }
        }
        //System.out.print(prev + "-" + number + "\n");

        int number_of_bins = 5;
        ArrayList<String[]>[] bins = new ArrayList[number_of_bins];
        for (bin = 0; bin < number_of_bins; bin++) {
            bins[bin] = new ArrayList();
            //System.out.println(bin + " " + bins[bin].size());
        }
        bin = 0;
        String id, suf;
        for (Node n : all.nodes.values()) {
            if (n.type.equals("protein")) {
                degree = n.downnbrs.size() + n.revnbrs.size() + n.upnbrs.size();
                if (degree > 0) {
                    try {
                        bin = get_bin(degree);
                        if (bin == number_of_bins) {
                            bin = number_of_bins;
                        }
                        String[] hgentry = {hugo_by_id.get(n.id.substring(1))[1], hugo_by_id.get(n.id.substring(1))[0]};
                        bins[bin].add(hgentry);
                        //all_hugo_ids.add(hugo_by_id.get(n.id.substring(1))[1]);
                    } catch (NullPointerException e) {

                    }
                }
            }

        }

        for (Node n : all.nodes.values()) {
            if (n.type.equals("miRNA")) {

                String[] hgentry = {n.id, n.id};
                all_mirnas_ids.add(hgentry);
            } else if (n.type.toLowerCase().equals("chemical") || n.type.toUpperCase().equals("SMALLMOLECULE")) {
                String[] hgentry = {n.id, n.id};
                all_chemicals_ids.add(hgentry);
            }
        }

        for (bin = 0; bin < number_of_bins; bin++) {
            //System.out.println(bin + " " + bins[bin].size());
        }

        String[] next_gene;
        String flag;
        bin = 0;
        BufferedWriter random_hg_file = new BufferedWriter(new FileWriter(folder + "RandomHitGenes" + "_" + thread_prefix));
        BufferedWriter random_hg_file_lists = new BufferedWriter(new FileWriter(folder + "RandomHitGenes_lists" + "_" + thread_prefix));
        List<Integer>[] already_encountered_p = new ArrayList[10];
        List<Integer> already_encountered_m = new ArrayList();
        for (bin = 0; bin < number_of_bins; bin++) {
            already_encountered_p[bin] = new ArrayList();
        }
        for (int j = count * thread; j > count * (thread - 1) + 1; j--) {
            random_hg_file_lists.write("HitList" + j + "\n");
            for (int i = count_list_size - 1; i >= 0; i--) {

                if ("protein".equals(all.nodes.get(hit_list.get(i)).type) || "gene".equals(all.nodes.get(hit_list.get(i)).type)) {
//                if (!hit_list.get(i).startsWith("MIR")) {
                    try {
                        int h = all.nodes.get(hit_list.get(i)).downnbrs.size();
                    } catch (NullPointerException e) {
                        System.out.println(hit_list.get(i));
                    }
                    degree = all.nodes.get(hit_list.get(i)).downnbrs.size() + all.nodes.get(hit_list.get(i)).upnbrs.size() + all.nodes.get(hit_list.get(i)).revnbrs.size();
                    bin = get_bin(degree);
                    if (bin == number_of_bins) {
                        bin = number_of_bins;
                    }
                    randomInt = randomGenerator.nextInt(bins[bin].size());
                    next_gene = bins[bin].get(randomInt);
                    flag = "p";
                } else if (all.nodes.get(hit_list.get(i)).type.equals("miRNA")) {
                    randomInt = randomGenerator.nextInt(all_mirnas_ids.size());
                    next_gene = all_mirnas_ids.get(randomInt);
                    flag = "m";
                } else if (all.nodes.get(hit_list.get(i)).type.toLowerCase().equals("chemical") || all.nodes.get(hit_list.get(i)).type.toUpperCase().equals("SMALLMOLECULE")) {
                    randomInt = randomGenerator.nextInt(all_chemicals_ids.size());
                    next_gene = all_chemicals_ids.get(randomInt);
                    flag = "c";
                } else {

                    System.out.println("Unknown hit gene type " + hit_list.get(i) + " type " + all.nodes.get(hit_list.get(i)).type);

                    continue;
                }
                if ((flag.equals("p") && already_encountered_p[bin].contains(randomInt)) || (flag.equals("m") && already_encountered_m.contains(randomInt))) {
                    random_hg_file_lists.write(next_gene[0] + "\t" + next_gene[1] + "\n");
                } else {
                    random_hg_file_lists.write(next_gene[0] + "\t" + next_gene[1] + "\n");
                    random_hg_file.write(next_gene[0] + "\t" + next_gene[1] + "\n");
                    if (flag.equals("p")) {
                        already_encountered_p[bin].add(randomInt);
                    } else {
                        already_encountered_m.add(randomInt);
                    }
                }
            }
        }
        random_hg_file.close();
        random_hg_file_lists.close();

        nutils.load_hitlist_and_finalimpl(all, folder + "RandomHitGenes" + "_" + thread_prefix, fplayers, dbutils.hugo);

//        System.out.println(prefix);
//        System.out.println(thread_prefix);
        nutils.find_pathway_for_list_BF_algorithm(all, dbutils.hugo, maxLength, folder + "random_paths_" + prefix + "_" + thread_prefix, prefix);
        putils.find_the_shortest_paths(folder + "random_paths_" + prefix + "_" + thread_prefix, "_shortest", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);

    }

    /**
     * Calculate p-values for paths for shuffled phenotype label
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     * @param random_paths_file File name for pathways from shuffled hit genes
     * @param ovrreps_paths_file File with paths
     * @param hitlists_file File with hit genes
     * @param hitlist_length Length of hit list
     * @param position_of_overreps_count Position of the centrality count in the
     * paths file
     * @param hugo_by_id HGNC ids map
     */
    public void calculate_p_values_phenotype_label_permutation(
            Network all,
            NetworkManager nutils,
            PathwayManager putils,
            CentralityManager outils,
            DBManager dbutils,
            String random_paths_file,
            String ovrreps_paths_file,
            String hitlists_file,
            int hitlist_length,
            int position_of_overreps_count,
            Map<String, String[]> hugo_by_id
    ) throws IOException {
        System.out.println("+++++++++++++calculate_p_values++++++++++++");
        BufferedReader rndm = new BufferedReader(new FileReader(random_paths_file));
        BufferedReader overrep = new BufferedReader(new FileReader(ovrreps_paths_file));
        BufferedReader hitlist = new BufferedReader(new FileReader(hitlists_file));
        BufferedWriter overrep_pvalue = new BufferedWriter(new FileWriter(ovrreps_paths_file + "_pvalues"));

        String s;
        String[] ss;
        List<String> tmp_list;
        String tmp_s;
        Map<String, List<String>> random_paths = new HashMap();
        Map<String, String> overreps = new HashMap();
        Map<String, Map<String, Integer>> overreps_counts = new HashMap();
        Map<String, List<String>> hitlists = new HashMap();
        Map<String, Integer> tmp_map;
        /*       while ((s = rndm.readLine()) != null) {
         ss = s.split(("\t"));
         if (random_paths.containsKey(ss[1])) {
         tmp_list = random_paths.get(ss[1]);
         tmp_list.add(s);
         } else {
         tmp_list = new ArrayList();
         tmp_list.add(s);
         }
         random_paths.put(ss[1], tmp_list);
         } */

        int j = 0;
        while ((s = overrep.readLine()) != null) {
            ss = s.split(("\t"));
            //System.out.println(ss[position_of_overreps_count]);
            if (Integer.parseInt(ss[position_of_overreps_count]) > 0 || (Integer.parseInt(ss[0])) >= 0 && Integer.parseInt(ss[position_of_overreps_count]) >= 0) {
                tmp_s = ss[1];
                for (int i = 2; i <= Integer.parseInt(ss[0]); i++) {
                    tmp_s = tmp_s + "\t" + ss[i];
                }
                overreps.put(tmp_s, s);
            }
        }

        int tmp_count;
        String hg_hugo_id;
        while ((s = rndm.readLine()) != null) {
            ss = s.split(("\t")); // ss[1] - hitgene
            try {
                if (ss[1].startsWith("hsa-") || ss[1].startsWith("-mir-")) {
                    hg_hugo_id = ss[1];
                } else if (ss[1].startsWith("CID:")) {
                    hg_hugo_id = ss[1];
                } else if (ss[1].startsWith("p")) {
                    hg_hugo_id = "p" + hugo_by_id.get(ss[1].substring(1))[1];
                    //System.out.println(hg_hugo_id);
                } else if (ss[1].startsWith("HGNC:")) {
                    hg_hugo_id = hugo_by_id.get(ss[1])[1];

                    //System.out.println(ss[1] + "\t" + hg_hugo_id);
                    //hg_hugo_id = "CAUTION";
                    //System.out.println("CAUTION");
                    //break;
                } else {
                    hg_hugo_id = ss[1];
//                    hg_hugo_id = "CAUTION";
//                    System.out.println(ss[1] + " CAUTION");
//                    continue;
                }
                for (String key : overreps.keySet()) {
                    if (s.contains(key)) {
                        if (overreps_counts.containsKey(key)) {
                            if (overreps_counts.get(key).containsKey(hg_hugo_id)) {
                                tmp_count = overreps_counts.get(key).get(hg_hugo_id) + 1;
                                overreps_counts.get(key).put(hg_hugo_id, tmp_count);
                            } else {
                                tmp_count = 1;
                                overreps_counts.get(key).put(hg_hugo_id, tmp_count);
                            }
                        } else {
                            tmp_map = new HashMap();
                            tmp_map.put(hg_hugo_id, 1);
                            overreps_counts.put(key, tmp_map);
                        }
                    }
                    /*else {
                     if (!overreps_counts.containsKey(key)) {
                     tmp_map = new HashMap();
                     tmp_map.put(hg_hugo_id, 0);
                     overreps_counts.put(key, tmp_map);
                     } else {
                     if (!overreps_counts.get(key).containsKey(hg_hugo_id)) {
                     tmp_map = overreps_counts.get(key);
                     tmp_map.put(hg_hugo_id, 0);
                     overreps_counts.put(key, tmp_map);
                     }
                     }
                     } */

                }
            } catch (NullPointerException e) {
                System.out.println(" --------- " + ss[1].substring(1));
            }
        }

        /*
         for (String key : overreps.keySet()) {
         tmp_map = new HashMap();
         for (String hg : random_paths.keySet()) {
         tmp_count = 0;
         for (String s1 : random_paths.get(hg)) {
         if (s1.contains(key)) {
         tmp_count++;
         }
         }
         try {
         tmp_map.put("p" + hugo_by_id.get(hg.substring(1))[1], tmp_count);
         } catch (NullPointerException e) {
         System.out.println(" --------- " + hg.substring(1));
         }
         }
         overreps_counts.put(key, tmp_map);
         j++;
         }
         */
        while ((s = hitlist.readLine()) != null) {
            if (s.contains("HitList")) {
                tmp_s = s;
                tmp_list = new ArrayList();
                for (int i = 0; i < hitlist_length; i++) {
                    s = hitlist.readLine();
                    tmp_list.add(s);
                }
                hitlists.put(tmp_s, tmp_list);
            }
        }
        double pvalue;
        int count_hitlists;
        int overreps_counts_value = 0;
        String id = "";
        for (String key : overreps_counts.keySet()) {
            overrep_pvalue.write(overreps.get(key) + "\t");
            count_hitlists = 0;
            for (String hl : hitlists.keySet()) {
                tmp_count = 0;
                for (String hg : hitlists.get(hl)) {
                    try {
                        id = returnGeneName(hg, hugo_by_id);
                        /*
                         if (hg.startsWith("MIR") && (hugo_by_id.get(hg)[2].startsWith("microRNA"))) {
                         if (hg.contains("LET")) {
                         id = "hsa-let-" + hg.substring(6).toLowerCase();
                         } else {
                         id = "hsa-mir-" + hg.substring(3).toLowerCase();
                         }
                         } else {
                         id = "p" + hg;// hit gene - gene                    
                         } */
                        //System.out.println("hg " + hg +" id " + id + " " + hl);
                        overreps_counts_value = overreps_counts.get(key).get(id);
                    } catch (NullPointerException e) {
                        overreps_counts_value = 0;
                    }
                    tmp_count = tmp_count + overreps_counts_value;
                }
                if (tmp_count >= Integer.parseInt(overreps.get(key).split("\t")[position_of_overreps_count])) {
                    count_hitlists++;
                }
            }
            pvalue = (double) count_hitlists / (double) hitlists.size();

            overrep_pvalue.write(pvalue
                    + "\n");
        }

        overrep_pvalue.close();
        rndm.close();
        overrep.close();
        hitlist.close();
    }

    public void adjust_pvalues(
            Network all,
            NetworkManager nutils,
            PathwayManager putils,
            CentralityManager outils,
            DBManager dbutils,
            String ovrreps_paths_file_pval,
            String ovrreps_paths_file_pval_adj,
            int min_centrality,
            int position_of_overreps_count,
            String mode
    ) throws IOException {
        //BufferedReader overrep_pvalue = new BufferedReader(new FileReader(ovrreps_paths_file + "_pvalues"));
        // BufferedWriter overrep_pvalueadj = new BufferedWriter(new FileWriter(ovrreps_paths_file + "_pvaluesadj"));
        BufferedReader overrep_pvalue = new BufferedReader(new FileReader(ovrreps_paths_file_pval));
        BufferedWriter overrep_pvalueadj = new BufferedWriter(new FileWriter(ovrreps_paths_file_pval_adj));
        String s;
        String[] ss;
        ArrayList<Double> pvals = new ArrayList();
        HashMap<Double, Double> pvalsadj = new HashMap();
        ArrayList<Integer> ranks = new ArrayList();

        while ((s = overrep_pvalue.readLine()) != null) {
            ss = s.split(("\t"));
            int centrality = 0;
            if (mode == "paths") {
                centrality = Integer.parseInt(ss[position_of_overreps_count]);

            }
            if (mode == "nodes") {
                centrality = Integer.parseInt(ss[2]);
            }
            if (centrality < min_centrality) {
                continue;
            }
            pvals.add(Double.parseDouble(ss[ss.length - 1]));
        }
        
        overrep_pvalue.close();
        Collections.sort(pvals);
        int rank = 0;
        for (int i = 0; i < pvals.size(); i++) {
            if (i == 0) {
                rank = 1;
                ranks.add(rank);
                continue;
            }
            if (pvals.get(i) == pvals.get(i - 1)) {
                ranks.add(rank);
            } else {
                rank += 1;
                ranks.add(rank);
            }
        }

        for (int i = 0; i < pvals.size(); i++) {
            Double pvaladj = pvals.get(i) / ranks.get(i) * pvals.size();

            pvalsadj.put(pvals.get(i), pvaladj);
        }

        int check = 1;

        while (true) {
            check = 0;
            for (int i = 0; i < pvals.size() - 1; i++) {
                if (pvalsadj.get(pvals.get(i)) > pvalsadj.get(pvals.get(i + 1))) {
                    pvalsadj.put(pvals.get(i), pvalsadj.get(pvals.get(i + 1)));
                    check = 1;
                }
            }
            if (check == 0) {
                    break;
            }

        }

        overrep_pvalue = new BufferedReader(new FileReader(ovrreps_paths_file_pval));
        while ((s = overrep_pvalue.readLine()) != null) {
            ss = s.split(("\t"));
            Double pval = Double.parseDouble(ss[ss.length - 1]);
            Double pvaladj = pvalsadj.get(pval);

            if (pvaladj > 1.0) {
                pvaladj = 1.0;
            }
            
            overrep_pvalueadj.write(s + "\t" + pvaladj + "\n");
        }
        overrep_pvalueadj.close();
        overrep_pvalue.close();
                
    }

    /**
     * Calculate p-values for nodes for shuffled phenotype label
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     * @param random_paths_file File name for pathways from shuffled hit genes
     * @param ovrreps_paths_file File with nodes
     * @param hitlists_file File with hit genes
     * @param hitlist_length Length of hit list
     * @param hugo_by_id HGNC ids map
     */
    public void calculate_p_values_phenotype_label_permutation_nodes(
            Network all,
            NetworkManager nutils,
            PathwayManager putils,
            CentralityManager outils,
            DBManager dbutils,
            String random_paths_file,
            String ovrreps_paths_file,
            String hitlists_file,
            int hitlist_length,
            Map<String, String[]> hugo_by_id
    ) throws IOException {
        System.out.println("+++++++++++++calculate_p_values_nodes++++++++++++");
        BufferedReader rndm = new BufferedReader(new FileReader(random_paths_file));
        BufferedReader overrep = new BufferedReader(new FileReader(ovrreps_paths_file));
        BufferedReader hitlist = new BufferedReader(new FileReader(hitlists_file));
        BufferedWriter overrep_pvalue = new BufferedWriter(new FileWriter(ovrreps_paths_file + "_pvalues"));

        String s;
        String[] ss;
        List<String> tmp_list;
        String tmp_s;
        Map<String, List<String>> random_paths = new HashMap();
        Map<String, String> hubs = new HashMap();
        Map<String, Map<String, Integer>> overreps_counts = new HashMap();
        Map<String, List<String>> hitlists = new HashMap();
        Map<String, Integer> tmp_map;
        int j = 0;
        while ((s = overrep.readLine()) != null) {
            ss = s.split(("\t"));
            if (Integer.parseInt(ss[2]) > 1) {
                hubs.put(ss[0], s);
            }
        }

        int tmp_count;
        String hg_hugo_id;
        boolean found;
        while ((s = rndm.readLine()) != null) {
            ss = s.split(("\t")); // ss[1] - hitgene
            try {
                if (ss[1].startsWith("hsa-") || ss[1].startsWith("-mir-")) {
                    hg_hugo_id = ss[1];
                } else if (ss[1].startsWith("CID:")) {
                    hg_hugo_id = ss[1];
                } else if (ss[1].startsWith("p")) {
                    hg_hugo_id = "p" + hugo_by_id.get(ss[1].substring(1))[1];
                    //System.out.println(hg_hugo_id);
                } else if (ss[1].startsWith("HGNC:")) {
                    hg_hugo_id = hugo_by_id.get(ss[1])[1];

                    //System.out.println(ss[1] + "\t" + hg_hugo_id);
                    //hg_hugo_id = "CAUTION";
                    //System.out.println("CAUTION");
                    //break;
                } else {
                    hg_hugo_id = ss[1];
//                    hg_hugo_id = "CAUTION";
//                    System.out.println(ss[1] + " CAUTION");
//                    continue;
                }

            } catch (NullPointerException e) {
                System.out.println(" --------- " + ss[1].substring(1));
                break;
            }
            for (String key : hubs.keySet()) {
                found = false;
                for (int k = 2; k < ss.length - 2; k++) {
                    try {
                        if (all.interactions.get(ss[k]).int1.id.equals(key)
                                || all.interactions.get(ss[k]).int2.id.equals(key)) {
                            found = true;
                            break;
                        }
                    } catch (NullPointerException e) {
                        System.out.println(" ++++ " + ss[k]);
                        break;
                    }
                }
                if (found) {
                    if (overreps_counts.containsKey(key)) {
                        if (overreps_counts.get(key).containsKey(hg_hugo_id)) {
                            tmp_count = overreps_counts.get(key).get(hg_hugo_id) + 1;
                            overreps_counts.get(key).put(hg_hugo_id, tmp_count);
                        } else {
                            tmp_count = 1;
                            overreps_counts.get(key).put(hg_hugo_id, tmp_count);
                        }
                    } else {
                        tmp_map = new HashMap();
                        tmp_map.put(hg_hugo_id, 1);
                        overreps_counts.put(key, tmp_map);
                    }
                }
            }

        }
        while ((s = hitlist.readLine()) != null) {
            if (s.contains("HitList")) {
                tmp_s = s;
                tmp_list = new ArrayList();
                for (int i = 0; i < hitlist_length; i++) {
                    s = hitlist.readLine();
                    tmp_list.add(s);
                }
                hitlists.put(tmp_s, tmp_list);
            }
        }
        double pvalue;
        int count_hitlists;
        String id = "";
        int overreps_counts_value = 0;
        for (String key : overreps_counts.keySet()) {
            overrep_pvalue.write(hubs.get(key) + "\t");
            count_hitlists = 0;
            for (String hl : hitlists.keySet()) {
                tmp_count = 0;
                for (String hg : hitlists.get(hl)) {
                    try {
                        id = returnGeneName(hg, hugo_by_id);
                        overreps_counts_value = overreps_counts.get(key).get(id);
                        //overreps_counts_value = overreps_counts.get(key).get("p" + hg);
                    } catch (NullPointerException e) {
                        overreps_counts_value = 0;
                    }
                    tmp_count = tmp_count + overreps_counts_value;
                }
                if (tmp_count >= Integer.parseInt(hubs.get(key).split("\t")[2])) {
                    count_hitlists++;
                }
            }
            pvalue = (double) count_hitlists / (double) hitlists.size();
            overrep_pvalue.write(pvalue + "\n");
        }

        overrep_pvalue.close();
        rndm.close();
        overrep.close();
        hitlist.close();
    }

    /**
     *
     * @param all
     * @param nutils
     * @param putils
     * @param outils
     * @param dbutils
     * @param random_paths_file
     * @param ovrreps_paths_file
     * @param hitlists_file
     * @param hitlist_length
     * @param position_of_overreps_count
     * @param hugo_by_id
     * @throws IOException
     */
    public void calculate_p_values_paths_random_networks(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, String random_paths_file, String ovrreps_paths_file, String hitlists_file, int hitlist_length, int position_of_overreps_count, Map<String, String[]> hugo_by_id) throws IOException {
        System.out.println("+++++++++++++calculate_p_values++++++++++++");
        BufferedReader rndm = new BufferedReader(new FileReader(random_paths_file));
        BufferedReader overrep = new BufferedReader(new FileReader(ovrreps_paths_file));
        BufferedReader hitlist = new BufferedReader(new FileReader(hitlists_file));
        BufferedWriter overrep_pvalue = new BufferedWriter(new FileWriter(ovrreps_paths_file + "_pvalues_2"));

        String s;
        String[] ss;
        List<String> tmp_list;
        String tmp_s;
        Map<String, List<String>> random_paths = new HashMap();
        Map<String, String> overreps = new HashMap();
        Map<String, Map<String, Integer>> overreps_counts = new HashMap();
        Map<String, List<String>> hitlists = new HashMap();
        Map<String, Integer> tmp_map;
        while ((s = rndm.readLine()) != null) {
            ss = s.split(("\t"));
            if (random_paths.containsKey(ss[1])) {
                tmp_list = random_paths.get(ss[1]);
                tmp_list.add(s);
            } else {
                tmp_list = new ArrayList();
                tmp_list.add(s);
            }
            random_paths.put(ss[1], tmp_list);
        }

        int j = 0;
        while ((s = overrep.readLine()) != null) {
            ss = s.split(("\t"));
            if (Integer.parseInt(ss[position_of_overreps_count]) > 3 || (Integer.parseInt(ss[0])) >= 3 && Integer.parseInt(ss[position_of_overreps_count]) > 2) {
                tmp_s = ss[1];
                for (int i = 2; i <= Integer.parseInt(ss[0]); i++) {
                    tmp_s = tmp_s + "\t" + ss[i];
                }
                overreps.put(tmp_s, s);
            }
        }

        int tmp_count;
        for (String key : overreps.keySet()) {
            tmp_map = new HashMap();
            for (String hg : random_paths.keySet()) {
                tmp_count = 0;
                for (String s1 : random_paths.get(hg)) {
                    if (s1.contains(key)) {
                        tmp_count++;
                    }
                }
                try {
                    tmp_map.put("p" + hugo_by_id.get(hg.substring(1))[1], tmp_count);
                } catch (NullPointerException e) {
                    System.out.println(" --------- " + hg.substring(1));
                }
            }
            overreps_counts.put(key, tmp_map);
            j++;
        }

        while ((s = hitlist.readLine()) != null) {
            if (s.contains("HitList")) {
                tmp_s = s;
                tmp_list = new ArrayList();
                for (int i = 0; i < hitlist_length; i++) {
                    s = hitlist.readLine();
                    tmp_list.add(s);
                }
                hitlists.put(tmp_s, tmp_list);
            }
        }
        double pvalue;
        int count_hitlists;
        for (String key : overreps_counts.keySet()) {
            overrep_pvalue.write(overreps.get(key) + "\t");
            count_hitlists = 0;
            for (String hl : hitlists.keySet()) {
                tmp_count = 0;
                for (String hg : hitlists.get(hl)) {
                    // System.out.println(overreps_counts.get(key).get("p" + hg));
                    try {
                        tmp_count = tmp_count + overreps_counts.get(key).get("p" + hg);
                    } catch (NullPointerException e) {

                    }
                }
                // System.out.println(tmp_count + "\t"+Integer.parseInt(overreps.get(key).split("\t")[position_of_overreps_count]));
                if (tmp_count >= Integer.parseInt(overreps.get(key).split("\t")[position_of_overreps_count])) {
                    count_hitlists++;
                }

            }
            pvalue = (double) count_hitlists / (double) hitlists.size();
            // System.out.println(pvalue);
            overrep_pvalue.write(pvalue + "\n");
        }
        overrep_pvalue.close();
        rndm.close();
        overrep.close();
        hitlist.close();
    }

    /**
     * Calculate p-values for paths for shuffled phenotype label for ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     * @param random_paths_file File name for pathways from shuffled hit genes
     * @param ovrreps_paths_file File with paths
     * @param hitlists_file File with hit genes
     * @param hitlist_length Length of hit list
     * @param position_of_overreps_count Position of the centrality count in the
     * paths file
     * @param hugo_by_id HGNC ids map
     */
    public void calculate_p_values_phenotype_label_permutation_ppi(
            Network all,
            NetworkManager nutils,
            PathwayManager putils,
            CentralityManager outils,
            DBManager dbutils,
            String random_paths_file,
            String ovrreps_paths_file,
            String hitlists_file,
            int hitlist_length,
            int position_of_overreps_count,
            Map<String, String[]> hugo_by_id
    ) throws IOException {
        System.out.println("+++++++++++++calculate_p_values++++++++++++");
        BufferedReader rndm = new BufferedReader(new FileReader(random_paths_file));
        BufferedReader overrep = new BufferedReader(new FileReader(ovrreps_paths_file));
        BufferedReader hitlist = new BufferedReader(new FileReader(hitlists_file));
        BufferedWriter overrep_pvalue = new BufferedWriter(new FileWriter(ovrreps_paths_file + "_pvalues"));

        String s;
        String[] ss;
        List<String> tmp_list;
        String tmp_s;
        Map<String, List<String>> random_paths = new HashMap();
        Map<String, String> overreps = new HashMap();
        Map<String, Map<String, Integer>> overreps_counts = new HashMap();
        Map<String, List<String>> hitlists = new HashMap();
        Map<String, Integer> tmp_map;
        //System.out.println("0");
        while ((s = rndm.readLine()) != null) {
            ss = s.split(("\t"));
            if (random_paths.containsKey(ss[1])) {
                tmp_list = random_paths.get(ss[1]);
                tmp_list.add(s);
            } else {
                tmp_list = new ArrayList();
                tmp_list.add(s);
            }
            random_paths.put(ss[1], tmp_list);
        }
        // System.out.println("1");
        int j = 0;
        while ((s = overrep.readLine()) != null) {
            ss = s.split(("\t"));
            if (Integer.parseInt(ss[position_of_overreps_count]) >= 0 || (Integer.parseInt(ss[0])) >= 0 && Integer.parseInt(ss[position_of_overreps_count]) >= 0) {
                tmp_s = ss[1];
                for (int i = 2; i <= Integer.parseInt(ss[0]); i++) {
                    tmp_s = tmp_s + "\t" + ss[i];
                }
                overreps.put(tmp_s, s);
                // System.out.println(tmp_s);                
            }
        }
        // int j=0;
        int tmp_count;
        // System.out.println("2");
        /*       for (String key : overreps.keySet()) {
         tmp_map = new HashMap();            
         for (String hg : random_paths.keySet()) {
         tmp_count = 0;
         for (String s1 : random_paths.get(hg)) {                    
         if (s1.contains(key)) {
         tmp_count++;
         // System.out.println( j + "\t"+s1);
         //j++;
         }
         if (s1.contains(CentralityManager.return_symmetric_path(key))) {
         tmp_count++;
         // System.out.println( j + "\t"+s1);
         //j++;
         }
         }                
         tmp_map.put("p" + hugo_by_id.get(hg.substring(1))[1], tmp_count);
         // System.out.println("p" + hugo_by_id.get(hg.substring(1))[1]);
         }
         overreps_counts.put(key, tmp_map);
         // System.out.println(j+"/" + overreps.keySet().size()+ "\t" +key);            
         j++;
         } */
        String hg_name;
        for (String key : overreps.keySet()) {
            overreps_counts.put(key, new HashMap());
        }
        for (String hg : random_paths.keySet()) {
            //System.out.println(j + "/" + random_paths.size());
            j++;
            //System.out.println(hg); hg.substring(1)
            hg_name = "p" + hugo_by_id.get(hg)[1];
            tmp_count = 0;
            for (String s1 : random_paths.get(hg)) {
                for (String key : overreps.keySet()) {
                    if (!overreps_counts.get(key).containsKey(hg_name)) {
                        overreps_counts.get(key).put(hg_name, 0);
                    }
                    if (s1.contains(key)) {
                        // System.out.println(overreps_counts.get(key).get(hg_name));
                        overreps_counts.get(key).put(hg_name, overreps_counts.get(key).get(hg_name) + 1);
                        // System.out.println(overreps_counts.get(key).get(hg_name));
                    }
                    if (s1.contains(CentralityManager.return_symmetric_path(key))) {
                        overreps_counts.get(key).put(hg_name, overreps_counts.get(key).get(hg_name) + 1);
                    }
                }
            }
        }

        //  System.out.println("3");
        while ((s = hitlist.readLine()) != null) {
            if (s.contains("HitList")) {
                tmp_s = s;
                tmp_list = new ArrayList();
                for (int i = 0; i < hitlist_length; i++) {
                    s = hitlist.readLine();
                    tmp_list.add(s);
                }
                hitlists.put(tmp_s, tmp_list);
            }
        }
        double pvalue;
        int count_hitlists;
        // System.out.println("4");
        for (String key : overreps_counts.keySet()) {
            overrep_pvalue.write(overreps.get(key) + "\t");
            count_hitlists = 0;
            for (String hl : hitlists.keySet()) {
                tmp_count = 0;
                for (String hg : hitlists.get(hl)) {
                    // System.out.println(overreps_counts.get(key).get("p" + hg));
                    try {
                        tmp_count = tmp_count + overreps_counts.get(key).get("p" + hg);
                    } catch (NullPointerException e) {

                    }
                }
                // System.out.println(tmp_count + "\t"+Integer.parseInt(overreps.get(key).split("\t")[position_of_overreps_count]));
                if (tmp_count >= Integer.parseInt(overreps.get(key).split("\t")[position_of_overreps_count])) {
                    count_hitlists++;
                }
            }
            pvalue = (double) count_hitlists / (double) hitlists.size();
            // System.out.println(pvalue);
            overrep_pvalue.write(pvalue + "\n");
        }
        overrep_pvalue.close();
        rndm.close();
        overrep.close();
        hitlist.close();
    }

    /**
     * Create random network with the same degree distribution
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public Network create_random_network(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils) {
        java.util.Random random_number = new java.util.Random();
        int randint = 0, randint2 = 0, step;
        String randid = "", randid2 = "", type, quality, dir, node_type;
        Network all_random = new Network();
        Interaction intr;
        Node node2;
        List<String> all_proteins = new ArrayList(), all_proteins_hard = new ArrayList();
        List<String> all_genes = new ArrayList(), all_genes_hard = new ArrayList();
        List<String> all_mirnas = new ArrayList(), all_mirnas_hard = new ArrayList();
        List<String> all_compounds = new ArrayList(), all_compounds_hard = new ArrayList();
        List<String> all_glycans = new ArrayList(), all_glycans_hard = new ArrayList();
        List<String> all_complexes = new ArrayList(), all_complexes_hard = new ArrayList();
        List<String> all_chemicals = new ArrayList(), all_chemicals_hard = new ArrayList();
        List<String> all_pfs = new ArrayList(), all_pfs_hard = new ArrayList();
        List<String> all_phts = new ArrayList(), all_phts_hard = new ArrayList();

        List<Node> sorted_direct = new ArrayList();
        List<Node> sorted_indirect = new ArrayList();

        //List<String> all_chemicals = new ArrayList();
        for (Node n : all.nodes.values()) {
            if (n.type.toLowerCase().equals("protein")) {
                all_proteins_hard.add(n.id);
            } else if (n.type.toLowerCase().equals("gene")) {
                all_genes_hard.add(n.id);
            } else if (n.type.toLowerCase().equals("compound")) {
                all_compounds_hard.add(n.id);
            } else if (n.type.toLowerCase().equals("glycan")) {
                all_glycans_hard.add(n.id);
            } else if (n.type.toLowerCase().equals("complex")) {
                all_complexes_hard.add(n.id);
            } else if (n.type.toLowerCase().equals("mirna")) {
                all_mirnas_hard.add(n.id);
            } else if (n.type.toLowerCase().equals("chemical") || n.type.toLowerCase().equals("smallmolecule")) {
                all_chemicals_hard.add(n.id);
            } else if (n.type.toLowerCase().equals("proteinfamily")) {
                all_pfs_hard.add(n.id);
            } else if (n.type.toLowerCase().equals("phenotype")) {
                all_phts_hard.add(n.id);
            } else {
                System.out.println(n.id + "\t" + n.type.toLowerCase());
            }
        }
        int index = 0;
        boolean found = false, add = false;
        Map<String, List<String>> down = new HashMap();
        Map<String, List<String>> up = new HashMap();
        Map<String, List<String>> rev = new HashMap();
        List<Interaction> not_found_interactions = new ArrayList();

        Node nn, nn2;
        boolean sort_list_add = false;
        int all_degree1, all_degree2, all_random_degree1, all_random_degree2;
        int sum_degree = 0, count_not_found = 0, sum_degree_random = 0;
        for (int j = 1; j > 0; j--) {
            all_random = new Network();

            all_proteins.addAll(all_proteins_hard);
            all_genes.addAll(all_genes_hard);
            all_mirnas.addAll(all_mirnas_hard);
            all_compounds.addAll(all_compounds_hard);
            all_glycans.addAll(all_glycans_hard);
            all_complexes.addAll(all_complexes_hard);
            for (String s : all.nodes.keySet()) {
                sort_list_add = false;
                nn = new Node(all.nodes.get(s).id, all.nodes.get(s).type, all.nodes.get(s).id_type, all.nodes.get(s).db_flag, all.nodes.get(s).ids);
                all_random.nodes.put(all.nodes.get(s).id, nn);
                sum_degree = sum_degree + all.nodes.get(s).downnbrs.size()
                        + all.nodes.get(s).upnbrs.size()
                        + all.nodes.get(s).revnbrs.size();

                if (!all.nodes.get(s).downnbrs.isEmpty()) {
                    if (sorted_direct.isEmpty()) {
                        sorted_direct.add(all.nodes.get(s));
                    }
                    for (int i = 0; i < sorted_direct.size(); i++) {
                        if (sorted_direct.get(i).downnbrs.size() <= all.nodes.get(s).downnbrs.size()) {
                            sorted_direct.add(i, all.nodes.get(s));
                            sort_list_add = true;
                            break;
                        }
                    }
                    if (!sort_list_add) {
                        sorted_direct.add(all.nodes.get(s));
                    }

                }
                sort_list_add = false;
                if (!all.nodes.get(s).revnbrs.isEmpty()) {
                    if (sorted_indirect.isEmpty()) {
                        sorted_indirect.add(all.nodes.get(s));
                    }
                    for (int i = 0; i < sorted_indirect.size(); i++) {
                        if (sorted_indirect.get(i).revnbrs.size() <= all.nodes.get(s).revnbrs.size()) {
                            sorted_indirect.add(i, all.nodes.get(s));
                            sort_list_add = true;
                            break;
                        }
                    }
                    if (!sort_list_add) {
                        sorted_indirect.add(all.nodes.get(s));
                    }

                }

            }
            System.out.println(" sorted direct size " + sorted_direct.size() + " sorted indirect size " + sorted_indirect.size() + " all nodes size " + all.nodes.size());

            for (Node n : sorted_direct) {
                for (Interaction it : n.downnbrs.values()) {
                    index++;
                    if (n.type.equals("gene") && it.int2.type.toLowerCase().equals("protein")) {
                        intr = new Interaction("R" + index, new ArrayList(), n, it.int2, it.type, "", new ArrayList(), it.quality, it.dir);
                        all_random.nodes.get(n.id).downnbrs.put(it.int2.id, intr);
                        all_random.nodes.get(it.int2.id).upnbrs.put(n.id, intr);
                        all_random.interactions.put("R" + index, intr);
                        continue;
                    }
                    type = it.type;
                    dir = it.dir;
                    quality = it.quality;
                    System.out.println(index);
                    step = 0;
                    add = false;
                    while (true) {
                        step++;
                        if (step < 10) {
                            randid = n.id;

                            if (it.int2.type.toLowerCase().equals("protein")) {
                                step = 0;
                                System.out.println("----------7");
                                while (true) {
                                    step++;
                                    if (step > 100) {
                                        add = false;
                                        count_not_found++;
                                        break;
                                    }

                                    if (all_proteins.size() == 0) {
                                        System.out.println("RRRRRRRR");
                                    }

                                    randint2 = random_number.nextInt(all_proteins.size());
                                    randid2 = all_proteins.get(randint2);
                                    if (randid.equals(randid2) || all_random.nodes.get(randid).downnbrs.containsKey(randid2)) {
                                        continue;
                                    }
                                    all_degree2 = all.nodes.get(randid2).upnbrs.size();
                                    all_random_degree2 = all_random.nodes.get(randid2).upnbrs.size();
                                    if (all.nodes.get(randid2).upnbrs.size() == all_random.nodes.get(randid2).upnbrs.size()) {
                                        all_proteins.remove(randid2);
                                    }
                                    if (all_random_degree2 < all_degree2) {
                                        add = true;
                                        break;
                                    }
                                }
                            }

                            if (it.int2.type.toLowerCase().equals("gene")) {
                                System.out.println("----------8");
                                while (true) {
                                    randint2 = random_number.nextInt(all_genes.size());
                                    randid2 = all_genes.get(randint2);
                                    if (randid.equals(randid2) || all_random.nodes.get(randid).downnbrs.containsKey(randid2)) {
                                        continue;
                                    }

                                    all_degree2 = all.nodes.get(randid2).upnbrs.size();
                                    all_random_degree2 = all_random.nodes.get(randid2).upnbrs.size();

                                    if (all.nodes.get(randid2).upnbrs.size() == all_random.nodes.get(randid2).upnbrs.size()) {
                                        all_genes.remove(randid2);
                                    }
                                    if (all_random_degree2 < all_degree2) {
                                        add = true;
                                        break;
                                    }
                                }
                            }

                            if (it.int2.type.toLowerCase().equals("mirna")) {
                                System.out.println("----------9");
                                while (true) {
                                    randint2 = random_number.nextInt(all_mirnas.size());
                                    randid2 = all_mirnas.get(randint2);
                                    if (randid.equals(randid2) || all_random.nodes.get(randid).downnbrs.containsKey(randid2)) {
                                        continue;
                                    }

                                    all_degree2 = all.nodes.get(randid2).upnbrs.size();
                                    all_random_degree2 = all_random.nodes.get(randid2).upnbrs.size();

                                    if (all.nodes.get(randid2).revnbrs.size() == all_random.nodes.get(randid2).revnbrs.size()
                                            && all.nodes.get(randid2).upnbrs.size() == all_random.nodes.get(randid2).upnbrs.size()
                                            && all.nodes.get(randid2).downnbrs.size() == all_random.nodes.get(randid2).downnbrs.size()) {
                                        all_mirnas.remove(randid2);
                                    }
                                    System.out.println(all_degree2 + " " + all_random_degree2 + " " + all_mirnas.size());

                                    if (all_random_degree2 < all_degree2) {
                                        add = true;
                                        break;
                                    }
                                }
                            }

                            if (it.int2.type.toLowerCase().equals("compound")) {
                                while (true) {
                                    randint2 = random_number.nextInt(all_compounds.size());
                                    randid2 = all_compounds.get(randint2);
                                    if (randid.equals(randid2) || all_random.nodes.get(randid).downnbrs.containsKey(randid2)) {
                                        continue;
                                    }
                                    System.out.println("----------10 " + randint2 + " " + randid2);

                                    all_degree2 = all.nodes.get(randid2).upnbrs.size();
                                    all_random_degree2 = all_random.nodes.get(randid2).upnbrs.size();

                                    System.out.println(all_degree2 + " " + all_random_degree2);
                                    if (all.nodes.get(randid2).revnbrs.size() == all_random.nodes.get(randid2).revnbrs.size()
                                            && all.nodes.get(randid2).upnbrs.size() == all_random.nodes.get(randid2).upnbrs.size()
                                            && all.nodes.get(randid2).downnbrs.size() == all_random.nodes.get(randid2).downnbrs.size()) {
                                        all_compounds.remove(randid2);
                                    }
                                    if (all_random_degree2 < all_degree2) {
                                        add = true;
                                        break;
                                    }
                                }
                            }

                            if (it.int2.type.toLowerCase().equals("glycan")) {
                                System.out.println("----------11");
                                while (true) {
                                    randint2 = random_number.nextInt(all_glycans.size());
                                    randid2 = all_glycans.get(randint2);
                                    if (randid.equals(randid2) || all_random.nodes.get(randid).downnbrs.containsKey(randid2)) {
                                        continue;
                                    }

                                    all_degree2 = all.nodes.get(randid2).upnbrs.size();
                                    all_random_degree2 = all_random.nodes.get(randid2).upnbrs.size();

                                    if (all.nodes.get(randid2).revnbrs.size() == all_random.nodes.get(randid2).revnbrs.size()
                                            && all.nodes.get(randid2).upnbrs.size() == all_random.nodes.get(randid2).upnbrs.size()
                                            && all.nodes.get(randid2).downnbrs.size() == all_random.nodes.get(randid2).downnbrs.size()) {
                                        all_glycans.remove(randid2);
                                    }
                                    if (all_random_degree2 < all_degree2) {
                                        add = true;
                                        break;
                                    }
                                }
                            }

                            if (it.int2.type.toLowerCase().equals("complex")) {
                                System.out.println("----------12");
                                while (true) {
                                    randint2 = random_number.nextInt(all_complexes.size());
                                    randid2 = all_complexes.get(randint2);
                                    if (randid.equals(randid2) || all_random.nodes.get(randid).downnbrs.containsKey(randid2)) {
                                        continue;
                                    }

                                    all_degree2 = all.nodes.get(randid2).upnbrs.size();
                                    all_random_degree2 = all_random.nodes.get(randid2).upnbrs.size();

                                    if (all.nodes.get(randid2).revnbrs.size() == all_random.nodes.get(randid2).revnbrs.size()
                                            && all.nodes.get(randid2).upnbrs.size() == all_random.nodes.get(randid2).upnbrs.size()
                                            && all.nodes.get(randid2).downnbrs.size() == all_random.nodes.get(randid2).downnbrs.size()) {
                                        all_complexes.remove(randid2);
                                    }
                                    if (all_random_degree2 < all_degree2) {
                                        add = true;
                                        break;
                                    }
                                }
                            }

                            if (it.int2.type.toLowerCase().equals("stimulus") || it.int2.type.toLowerCase().equals("chemical")
                                    || it.int2.type.toLowerCase().equals("smallmolecule")
                                    || it.int2.type.toLowerCase().equals("phenotype")
                                    || it.int2.type.toLowerCase().equals("proteinfamily")) {
                                continue;
                            }

                            if (!all_random.nodes.get(randid).downnbrs.containsKey(randid2)) {
                                add = true;
                                break;
                            }

                        }

                        if (add || (step >= 20)) {
                            if (step >= 20) {
                                // count_not_found++;
                            }
                            break;
                        }

                    }
                    if (add) {
                        intr = new Interaction("R" + index, new ArrayList(), all.nodes.get(randid), all.nodes.get(randid2), type, "", new ArrayList(), quality, dir);

                        if (all_random.nodes.get(randid).downnbrs.containsKey(randid2) || all_random.nodes.get(randid2).upnbrs.containsKey(randid)) {
                            System.out.println("AAAAAAACCCCCCCCCCHHHHHHHHHHHTTTTTTTTTTTUUUUUUUUUNNNNNNNNNNNNGGGGGG2");
                            break;
                        }

                        all_random.nodes.get(randid).downnbrs.put(randid2, intr);
                        all_random.nodes.get(randid2).upnbrs.put(randid, intr);
                        all_random.interactions.put("R" + index, intr);
                    }
                }

            }

            all_proteins = new ArrayList();
            all_genes = new ArrayList();
            all_mirnas = new ArrayList();
            all_compounds = new ArrayList();
            all_glycans = new ArrayList();
            all_complexes = new ArrayList();
            all_proteins.addAll(all_proteins_hard);
            all_genes.addAll(all_genes_hard);
            all_mirnas.addAll(all_mirnas_hard);
            all_compounds.addAll(all_compounds_hard);
            all_glycans.addAll(all_glycans_hard);
            all_complexes.addAll(all_complexes_hard);

            for (Node n : sorted_indirect) {
                for (Interaction it : n.revnbrs.values()) {
                    if (it.int1.id.equals(n.id)) {
                        node2 = it.int2;
                    } else {
                        node2 = it.int1;
                    }
                    if (sorted_indirect.indexOf(node2) > sorted_indirect.indexOf(n)) {
                        continue;
                    }
                    index++;
                    type = it.type;
                    dir = it.dir;
                    quality = it.quality;
                    //   System.out.println(index);
                    step = 0;
                    add = false;
                    while (true) {
                        step++;
                        if (step < 10) {
                            randid = n.id;

                            if (node2.type.toLowerCase().equals("protein")) {
                                step = 0;
                                //      System.out.println("----------107");
                                while (true) {
                                    step++;
                                    if (step > 100) {
                                        add = false;
                                        count_not_found++;
                                        break;
                                    }

                                    if (all_proteins.size() == 0) {
                                        System.out.println("RRRRRRRR");
                                    }

                                    randint2 = random_number.nextInt(all_proteins.size());
                                    randid2 = all_proteins.get(randint2);
                                    if (randid.equals(randid2) || all_random.nodes.get(randid).revnbrs.containsKey(randid2)) {
                                        continue;
                                    }
                                    all_degree2 = all.nodes.get(randid2).revnbrs.size();
                                    all_random_degree2 = all_random.nodes.get(randid2).revnbrs.size();
                                    if (all.nodes.get(randid2).revnbrs.size() == all_random.nodes.get(randid2).revnbrs.size()) {
                                        all_proteins.remove(randid2);
                                    }
                                    if (all_random_degree2 < all_degree2) {
                                        add = true;
                                        break;
                                    }
                                }
                            }

                            if (node2.type.toLowerCase().equals("gene")) {
                                System.out.println("----------108");
                                while (true) {
                                    randint2 = random_number.nextInt(all_genes.size());
                                    randid2 = all_genes.get(randint2);
                                    if (randid.equals(randid2) || all_random.nodes.get(randid).revnbrs.containsKey(randid2)) {
                                        continue;
                                    }

                                    all_degree2 = all.nodes.get(randid2).revnbrs.size();
                                    all_random_degree2 = all_random.nodes.get(randid2).revnbrs.size();

                                    if (all.nodes.get(randid2).revnbrs.size() == all_random.nodes.get(randid2).revnbrs.size()) {
                                        all_genes.remove(randid2);
                                    }
                                    if (all_random_degree2 < all_degree2) {
                                        add = true;
                                        break;
                                    }
                                }
                            }

                            if (node2.type.toLowerCase().equals("mirna")) {
                                System.out.println("----------109");
                                while (true) {
                                    randint2 = random_number.nextInt(all_mirnas.size());
                                    randid2 = all_mirnas.get(randint2);
                                    if (randid.equals(randid2) || all_random.nodes.get(randid).revnbrs.containsKey(randid2)) {
                                        continue;
                                    }

                                    all_degree2 = all.nodes.get(randid2).revnbrs.size();
                                    all_random_degree2 = all_random.nodes.get(randid2).revnbrs.size();

                                    if (all.nodes.get(randid2).revnbrs.size() == all_random.nodes.get(randid2).revnbrs.size()) {
                                        all_mirnas.remove(randid2);
                                    }

                                    if (all_random_degree2 < all_degree2) {
                                        add = true;
                                        break;
                                    }
                                }
                            }

                            if (node2.type.toLowerCase().equals("compound")) {
                                System.out.println("----------110");

                                while (true) {
                                    randint2 = random_number.nextInt(all_compounds.size());
                                    randid2 = all_compounds.get(randint2);
                                    if (randid.equals(randid2) || all_random.nodes.get(randid).revnbrs.containsKey(randid2)) {
                                        continue;
                                    }

                                    all_degree2 = all.nodes.get(randid2).revnbrs.size();
                                    all_random_degree2 = all_random.nodes.get(randid2).revnbrs.size();

                                    System.out.println(all_degree2 + " " + all_random_degree2);
                                    if (all.nodes.get(randid2).revnbrs.size() == all_random.nodes.get(randid2).revnbrs.size()) {
                                        all_compounds.remove(randid2);
                                    }
                                    if (all_random_degree2 < all_degree2) {
                                        add = true;
                                        break;
                                    }
                                }
                            }

                            if (node2.type.toLowerCase().equals("glycan")) {
                                System.out.println("----------111");
                                while (true) {
                                    randint2 = random_number.nextInt(all_glycans.size());
                                    randid2 = all_glycans.get(randint2);
                                    if (randid.equals(randid2) || all_random.nodes.get(randid).revnbrs.containsKey(randid2)) {
                                        continue;
                                    }

                                    all_degree2 = all.nodes.get(randid2).revnbrs.size();
                                    all_random_degree2 = all_random.nodes.get(randid2).revnbrs.size();

                                    if (all.nodes.get(randid2).revnbrs.size() == all_random.nodes.get(randid2).revnbrs.size()) {
                                        all_glycans.remove(randid2);
                                    }
                                    if (all_random_degree2 < all_degree2) {
                                        add = true;
                                        break;
                                    }
                                }
                            }

                            if (node2.type.toLowerCase().equals("complex")) {
                                System.out.println("----------112");
                                while (true) {
                                    randint2 = random_number.nextInt(all_complexes.size());
                                    randid2 = all_complexes.get(randint2);
                                    if (randid.equals(randid2) || all_random.nodes.get(randid).revnbrs.containsKey(randid2)) {
                                        continue;
                                    }

                                    all_degree2 = all.nodes.get(randid2).revnbrs.size();
                                    all_random_degree2 = all_random.nodes.get(randid2).revnbrs.size();

                                    if (all.nodes.get(randid2).revnbrs.size() == all_random.nodes.get(randid2).revnbrs.size()) {
                                        all_complexes.remove(randid2);
                                    }
                                    if (all_random_degree2 < all_degree2) {
                                        add = true;
                                        break;
                                    }
                                }
                            }

                            if (node2.type.toLowerCase().equals("stimulus") || it.int2.type.toLowerCase().equals("chemical")
                                    || it.int2.type.toLowerCase().equals("smallmolecule")
                                    || it.int2.type.toLowerCase().equals("phenotype")
                                    || it.int2.type.toLowerCase().equals("proteinfamily")) {
                                continue;
                            }

                            if (!all_random.nodes.get(randid).revnbrs.containsKey(randid2)) {
                                add = true;
                                break;
                            }

                        }

                        if (add || (step >= 20)) {
                            if (step >= 20) {
                                // count_not_found++;
                            }
                            break;
                        }

                    }
                    if (add) {
                        intr = new Interaction("R" + index, new ArrayList(), all.nodes.get(randid), all.nodes.get(randid2), type, "", new ArrayList(), quality, dir);

                        if (all_random.nodes.get(randid).revnbrs.containsKey(randid2) || all_random.nodes.get(randid2).revnbrs.containsKey(randid)) {
                            System.out.println("AAAAAAACCCCCCCCCCHHHHHHHHHHHTTTTTTTTTTTUUUUUUUUUNNNNNNNNNNNNGGGGGG2");
                            break;
                        }

                        all_random.nodes.get(randid).revnbrs.put(randid2, intr);
                        all_random.nodes.get(randid2).revnbrs.put(randid, intr);
                        all_random.interactions.put("R" + index, intr);
                    }
                }

            }

            for (String w : all_random.nodes.keySet()) {
                sum_degree_random = sum_degree_random + all_random.nodes.get(w).downnbrs.size()
                        + all_random.nodes.get(w).upnbrs.size()
                        + all_random.nodes.get(w).revnbrs.size();
            }

            sum_degree = sum_degree / 2;
            sum_degree_random = sum_degree_random / 2;
            System.out.println("interactions all " + all.interactions.size() + " random " + all_random.interactions.size());
            System.out.println("sum degree all " + sum_degree + " sum degree random " + sum_degree_random);
            System.out.println("count_not_found " + count_not_found);
            for (String s : all.nodes.keySet()) {
                if (all.nodes.get(s).downnbrs.size() != all_random.nodes.get(s).downnbrs.size()
                        || all.nodes.get(s).upnbrs.size() != all_random.nodes.get(s).upnbrs.size() //|| all.nodes.get(s).revnbrs.size() != all_random.nodes.get(s).revnbrs.size()
                        ) {
                    System.out.println(s
                            + " " + all.nodes.get(s).downnbrs.size() + " " + all_random.nodes.get(s).downnbrs.size()
                            + " " + all.nodes.get(s).upnbrs.size() + " " + all_random.nodes.get(s).upnbrs.size()
                            + " " + all.nodes.get(s).revnbrs.size() + " " + all_random.nodes.get(s).revnbrs.size());
                    for (String t : all.nodes.get(s).downnbrs.keySet()) {
                        if (!all_random.nodes.get(s).downnbrs.containsKey(t)) {
                            // System.out.println(s + "-" +t); 
                        }

                    }

                }
                //System.out.println(s + "\t" + all.nodes.get(s).downnbrs.size() + "\t" + all_random.nodes.get(s).downnbrs.size());
            }
        }
        return all_random;
    }

    /**
     * Calculate p-values for nodes on random network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     * @param file_name_random Name for file to store information about pathways
     * built on the random network
     * @param count Number of shuffling steps
     * @param hubs_file File with nodes and centrality scores
     * @param overreps_file File with paths and centrality scores
     */
    public void calculate_p_values_nodes_random_networks(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, String file_name_random, int count, String hubs_file, String overreps_file) throws IOException {
        BufferedReader rd = new BufferedReader(new FileReader(hubs_file));
        BufferedWriter wr = new BufferedWriter(new FileWriter(hubs_file + "_pvalues"));
        Map<String, String[]> hubs = new HashMap();
        Map<String, Integer> hubs_count = new HashMap();
        int tmp_int;
        String s;
        String[] ss;
        float pvalue;

        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            hubs.put(ss[0], ss);
            hubs_count.put(ss[0], 0);
        }
        rd.close();
        for (int j = count; j > 0; j--) {
            Network all_random = create_random_network(all, nutils, putils, outils, dbutils);
            nutils.find_pathway_for_list_BF_algorithm(all_random, dbutils.hugo, 6, file_name_random, "RD");
            putils.find_the_shortest_paths(file_name_random, "_top_ranked", all_random, nutils.hg, nutils.fpl, 0, 0, 0, 0);
            putils.add_connectivity_to_pathways(all_random, file_name_random + "_top_ranked", file_name_random + "_top_ranked_with_connectivity");
            outils.calculate_centrality_scores_for_paths(all, file_name_random + "_top_ranked", file_name_random + "_top_ranked__unique_overrepr", 1, 4);
            outils.add_hgnc_symbols_to_paths(file_name_random + "_top_ranked__unique_overrepr", file_name_random + "_top_ranked__unique_overrepr_w_names", all_random, dbutils.hugo_by_id);
            outils.calculate_paths_degree(file_name_random + "_top_ranked__unique_overrepr_w_names", file_name_random + "_top_ranked_with_connectivity", file_name_random + "_top_ranked__unique_overrepr_w_names_pthconn");
            outils.calculate_centrality_scores_for_nodes(all_random, dbutils.hugo_by_id, nutils, file_name_random + "_top_ranked", file_name_random + "_top_ranked_hubs");
            rd = new BufferedReader(new FileReader(file_name_random + "_top_ranked_hubs"));

            while ((s = rd.readLine()) != null) {
                ss = s.split("\t");
                if (hubs.containsKey(ss[0]) && Integer.parseInt(hubs.get(ss[0])[2]) > Integer.parseInt(ss[2])) {
                    tmp_int = hubs_count.get(ss[0]);
                    hubs_count.put(ss[0], tmp_int + 1);
                }

            }
            rd.close();

        }

        for (String t : hubs.keySet()) {
            pvalue = (float) hubs_count.get(t) / count;
            wr.write(hubs.get(t)[0] + "\t" + hubs.get(t)[1] + "\t" + hubs.get(t)[2] + "\t" + pvalue + "\n");
        }

        wr.close();
    }

    /**
     * Calculate p-values on for nodes random ppi network
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     * @param file_name_random File name to store information about pathways
     * built on the random network
     * @param count Number of shuffling steps
     * @param hubs_file File with nodes and centrality scores
     * @param overreps_file File with paths and centrality scores
     */
    public void calculate_p_values_nodes_random_networks_ppi(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, String file_name_random, int count, String hubs_file, String overreps_file) throws IOException {
        BufferedReader rd = new BufferedReader(new FileReader(hubs_file));
        BufferedWriter wr = new BufferedWriter(new FileWriter(hubs_file + "_pvalues"));
        Map<String, String[]> hubs = new HashMap();
        Map<String, Integer> hubs_count = new HashMap();
        int tmp_int;
        String s;
        String[] ss;
        float pvalue;

        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            hubs.put(ss[0], ss);
            hubs_count.put(ss[0], 0);
        }
        rd.close();
        for (int j = count; j > 0; j--) {
            Network all_random = create_random_network(all, nutils, putils, outils, dbutils);
            nutils.find_pathway_for_list_BF_algorithm_ppi(all_random, dbutils.hugo, 6, file_name_random, "RD");
            putils.find_the_shortest_paths(file_name_random, "_top_ranked", all_random, nutils.hg, nutils.fpl, 0, 0, 0, 0);
            putils.add_connectivity_to_pathways(all_random, file_name_random + "_top_ranked", file_name_random + "_top_ranked_with_connectivity");
            outils.calculate_centrality_scores_for_paths_ppi(all_random, file_name_random + "_top_ranked", file_name_random + "_top_ranked__unique_overrepr", 1, 4);
            outils.add_hgnc_symbols_to_paths(file_name_random + "_top_ranked__unique_overrepr", file_name_random + "_top_ranked__unique_overrepr_w_names", all_random, dbutils.hugo_by_id);
            outils.calculate_paths_degree(file_name_random + "_top_ranked__unique_overrepr_w_names", file_name_random + "_top_ranked_with_connectivity", file_name_random + "_top_ranked__unique_overrepr_w_names_pthconn");
            outils.calculate_centrality_scores_for_nodes(all_random, dbutils.hugo_by_id, nutils, file_name_random + "_top_ranked", file_name_random + "_top_ranked_hubs");

            rd = new BufferedReader(new FileReader(file_name_random + "_top_ranked_hubs"));
            while ((s = rd.readLine()) != null) {
                ss = s.split("\t");
                if (hubs.containsKey(ss[0]) && (Integer.parseInt(hubs.get(ss[0])[2]) <= Integer.parseInt(ss[2]))) {
                    tmp_int = hubs_count.get(ss[0]);
                    hubs_count.put(ss[0], tmp_int + 1);
                }

            }
            rd.close();
        }

        for (String t : hubs.keySet()) {
            pvalue = (float) hubs_count.get(t) / count;
            wr.write(hubs.get(t)[0] + "\t" + hubs.get(t)[1] + "\t" + hubs.get(t)[2] + "\t" + pvalue + "\n");
        }

        wr.close();
    }

    private void check_distributions(Network all, Map<String, String[]> hugo_by_id) {
        int number_of_bins = 5;
        List<String> all_hugo_ids = new ArrayList();
        List<String> all_mirnas_ids = new ArrayList();
        ArrayList<String>[] bins = new ArrayList[number_of_bins];
        int degree = 0;
        int bin = 0;
        String id = "", id2 = "", suf = "";
        int mirna_count = 0, mirna_con_count = 0, mirna_hugo_count = 0;
        for (Node n : all.nodes.values()) {
            if (n.type.equals("protein")) {
                degree = n.downnbrs.size() + n.revnbrs.size() + n.upnbrs.size();
                if (degree > 0) {
                    try {
                        bin = get_bin(degree);
                        if (bin == number_of_bins) {
                            bin = number_of_bins;
                        }
                        bins[bin].add(hugo_by_id.get(n.id.substring(1))[1]);
                        //all_hugo_ids.add(hugo_by_id.get(n.id.substring(1))[1]);
                    } catch (NullPointerException e) {

                    }
                }
            }
            if (n.type.toLowerCase().equals("miRNA".toLowerCase())) {
                mirna_count++;
                degree = n.downnbrs.size() + n.revnbrs.size() + n.upnbrs.size();
                if (degree > 0) {
                    mirna_con_count++;
                }
                id = n.id.replace("-3p", "").replace("-5p", "").replace("*", "");
                if (id.charAt(id.length() - 2) == '-' && Character.isDigit(id.charAt(id.length() - 1)) && Character.isDigit(id.charAt(id.length() - 3))) {
                    suf = "-" + id.substring(id.length() - 1);
                    id = id.substring(0, id.length() - 2);
                } else {
                    suf = "";
                }
                if (id.contains("let-")) {
                    id2 = "MIR" + id.replace("-", "").substring(3).toUpperCase() + suf;
                    all_mirnas_ids.add(id2);
                    //System.out.println(id + " " + id2);
                } else {
                    id2 = id.replace("-", "").substring(3).toUpperCase() + suf;
                    all_mirnas_ids.add(id2);
                    //System.out.println(id + " " + id2);
                }

            }
        }

        for (String t : hugo_by_id.keySet()) {

            if (hugo_by_id.get(t)[1].startsWith("MIR")) {
                //System.out.println(hugo_by_id.get(t)[1]);
                mirna_hugo_count++;
            }
        }
        System.out.println(mirna_hugo_count + " mirna_hugo_count " + mirna_count + " mirna_count " + mirna_con_count + " mirna_con_count");
    }

    private int get_bin(int degree) {
        boolean found = false;
        int i;
        int[] bins = new int[]{1, 3, 5, 10, 100};
        for (i = 0; i <= bins.length - 2; i++) {
            if (degree >= bins[i] && degree < bins[i + 1]) {
                return i;
            }
        }
        return i;

    }

    private String returnGeneName(String s, Map<String, String[]> hugo_by_id) {
        String index = "";
        String gene;
        if (s.startsWith("MIR")) {
            if (s.contains("LET")) {
                index = s.substring(6).toLowerCase();
                if (index.length() >= 2 && Character.isDigit(index.charAt(index.length() - 1)) && !Character.isDigit(index.charAt(index.length() - 2))) {
                    index = index.substring(0, index.length() - 2) + "-" + index.charAt(index.length() - 1);
                }
                gene = "hsa-let-" + index;
            } else {
                index = s.substring(3).toLowerCase();
                if (index.length() >= 2 && Character.isDigit(index.charAt(index.length() - 1)) && !Character.isDigit(index.charAt(index.length() - 2))) {
                    index = index.substring(0, index.length() - 2) + "-" + index.charAt(index.length() - 1);
                }
                gene = "hsa-mir-" + index;

            }
        } else {
            gene = "p" + hugo_by_id.get(s)[0];// hit gene - gene
        }
        return gene;
    }

}
/*


 if (it.int2.type.toLowerCase().equals("protein")) {
 while (true) {
 randint = random_number.nextInt(all_proteins.size());
 randid = all_proteins.get(randint);
 try {
 if (up.get(s).contains(randid) || up.get(randid).contains(s)) {
 continue;
 } else {
 break;
 }
 } catch (NullPointerException e) {
 break;
 }
 }
 }
 if (it.int2.type.toLowerCase().equals("gene")) {
 while (true) {
 randint = random_number.nextInt(all_genes.size());
 randid = all_genes.get(randint);
 try {
 if (up.get(s).contains(randid) || up.get(randid).contains(s)) {
 continue;
 } else {
 break;
 }
 } catch (NullPointerException e) {
 break;
 }
 }
 }
 if (it.int2.type.toLowerCase().equals("mirna")) {
 while (true) {
 randint = random_number.nextInt(all_mirnas.size());
 randid = all_mirnas.get(randint);
 try {
 if (up.get(s).contains(randid) || up.get(randid).contains(s)) {
 continue;
 } else {
 break;
 }
 } catch (NullPointerException e) {
 break;
 }
 }
 }
 if (it.int2.type.toLowerCase().equals("compound")) {
 while (true) {
 randint = random_number.nextInt(all_compounds.size());
 randid = all_compounds.get(randint);
 try {
 if (up.get(s).contains(randid) || up.get(randid).contains(s)) {
 continue;
 } else {
 break;
 }
 } catch (NullPointerException e) {
 break;
 }
 }
 }
 if (it.int2.type.toLowerCase().equals("glycan")) {
 while (true) {
 randint = random_number.nextInt(all_glycans.size());
 randid = all_glycans.get(randint);
 try {
 if (up.get(s).contains(randid) || up.get(randid).contains(s)) {
 continue;
 } else {
 break;
 }
 } catch (NullPointerException e) {
 break;
 }
 }
 }
 if (it.int2.type.toLowerCase().equals("complex")) {
 while (true) {
 randint = random_number.nextInt(all_complexes.size());
 randid = all_complexes.get(randint);
 try {
 if (up.get(s).contains(randid) || up.get(randid).contains(s)) {
 continue;
 } else {
 break;
 }
 } catch (NullPointerException e) {
 break;
 }
 }
 }
 if (it.int2.type.toLowerCase().equals("stimulus") || it.int2.type.toLowerCase().equals("chemical")
 || it.int2.type.toLowerCase().equals("smallmolucule")
 || it.int2.type.toLowerCase().equals("phenotype")
 || it.int2.type.toLowerCase().equals("proteinfamily")) {
 continue;
 }





 if (it.int2.type.toLowerCase().equals("protein")) {
 while (true) {
 randint = random_number.nextInt(all_proteins.size());
 randid = all_proteins.get(randint);
 try {
 if (down.get(s).contains(randid) || down.get(randid).contains(s)) {
 continue;
 } else {
 break;
 }
 } catch (NullPointerException e) {
 break;
 }
 }
 }
 if (it.int2.type.toLowerCase().equals("gene")) {
 while (true) {
 randint = random_number.nextInt(all_genes.size());
 randid = all_genes.get(randint);
 try {
 if (down.get(s).contains(randid) || down.get(randid).contains(s)) {
 continue;
 } else {
 break;
 }
 } catch (NullPointerException e) {
 break;
 }
 }
 }
 if (it.int2.type.toLowerCase().equals("mirna")) {
 while (true) {
 randint = random_number.nextInt(all_mirnas.size());
 randid = all_mirnas.get(randint);
 try {
 if (down.get(s).contains(randid) || down.get(randid).contains(s)) {
 continue;
 } else {
 break;
 }
 } catch (NullPointerException e) {
 break;
 }
 }
 }
 if (it.int2.type.toLowerCase().equals("compound")) {
 while (true) {
 randint = random_number.nextInt(all_compounds.size());
 randid = all_compounds.get(randint);
 try {
 if (down.get(s).contains(randid) || down.get(randid).contains(s)) {
 continue;
 } else {
 break;
 }
 } catch (NullPointerException e) {
 break;
 }
 }
 }
 if (it.int2.type.toLowerCase().equals("glycan")) {
 while (true) {
 randint = random_number.nextInt(all_glycans.size());
 randid = all_glycans.get(randint);
 try {
 if (down.get(s).contains(randid) || down.get(randid).contains(s)) {
 continue;
 } else {
 break;
 }
 } catch (NullPointerException e) {
 break;
 }
 }
 }
 if (it.int2.type.toLowerCase().equals("complex")) {
 while (true) {
 randint = random_number.nextInt(all_complexes.size());
 randid = all_complexes.get(randint);
 try {
 if (down.get(s).contains(randid) || down.get(randid).contains(s)) {
 continue;
 } else {
 break;
 }
 } catch (NullPointerException e) {
 break;
 }
 }
 }
 if (it.int2.type.toLowerCase().equals("stimulus") || it.int2.type.toLowerCase().equals("chemical")
 || it.int2.type.toLowerCase().equals("smallmolucule")
 || it.int2.type.toLowerCase().equals("phenotype")
 || it.int2.type.toLowerCase().equals("proteinfamily")) {
 continue;
 }


 */
 /*   if (all.nodes.get(randid2).revnbrs.size() == all_random.nodes.get(randid2).revnbrs.size()
 && all.nodes.get(randid2).upnbrs.size() == all_random.nodes.get(randid2).upnbrs.size()
 && all.nodes.get(randid2).downnbrs.size() == all_random.nodes.get(randid2).downnbrs.size()) {
 all_mirnas.remove(randid2);
 } */

 /*



 if (it.int1.type.toLowerCase().equals("protein")) {
 System.out.println("----------1");
 while (true) {
 randint = random_number.nextInt(all_proteins.size());
 randid = all_proteins.get(randint);
 if (dir.equals("indirect") || type.equals("kegg_reaction_reversible")) {
 all_degree1 = all.nodes.get(randid).revnbrs.size();
 all_random_degree1 = all_random.nodes.get(randid).revnbrs.size();
 } else {
 all_degree1 = all.nodes.get(randid).downnbrs.size();
 all_random_degree1 = all_random.nodes.get(randid).downnbrs.size();
 }
 if (all.nodes.get(randid).revnbrs.size() == all_random.nodes.get(randid).revnbrs.size()
 && all.nodes.get(randid).upnbrs.size() == all_random.nodes.get(randid).upnbrs.size()
 && all.nodes.get(randid).downnbrs.size() == all_random.nodes.get(randid).downnbrs.size()) {
 all_proteins.remove(randid);
 }
 if (all_random_degree1 < all_degree1) {
 break;
 }
 }
 }

 if (it.int1.type.toLowerCase().equals("gene")) {
 System.out.println("----------2");
 while (true) {
 randint = random_number.nextInt(all_genes.size());
 randid = all_genes.get(randint);
 if (dir.equals("indirect") || type.equals("kegg_reaction_reversible")) {
 all_degree1 = all.nodes.get(randid).revnbrs.size();
 all_random_degree1 = all_random.nodes.get(randid).revnbrs.size();
 } else {
 all_degree1 = all.nodes.get(randid).downnbrs.size();
 all_random_degree1 = all_random.nodes.get(randid).downnbrs.size();
 }
 if (all.nodes.get(randid).revnbrs.size() == all_random.nodes.get(randid).revnbrs.size()
 && all.nodes.get(randid).upnbrs.size() == all_random.nodes.get(randid).upnbrs.size()
 && all.nodes.get(randid).downnbrs.size() == all_random.nodes.get(randid).downnbrs.size()) {
 all_genes.remove(randid);
 }
 if (all_random_degree1 < all_degree1) {
 break;
 }
 }
 }

 if (it.int1.type.toLowerCase().equals("mirna")) {
 System.out.println("----------3");
 while (true) {
 randint = random_number.nextInt(all_mirnas.size());
 randid = all_mirnas.get(randint);
 if (dir.equals("indirect") || type.equals("kegg_reaction_reversible")) {
 all_degree1 = all.nodes.get(randid).revnbrs.size();
 all_random_degree1 = all_random.nodes.get(randid).revnbrs.size();
 } else {
 all_degree1 = all.nodes.get(randid).downnbrs.size();
 all_random_degree1 = all_random.nodes.get(randid).downnbrs.size();
 }
 if (all.nodes.get(randid).revnbrs.size() == all_random.nodes.get(randid).revnbrs.size()
 && all.nodes.get(randid).upnbrs.size() == all_random.nodes.get(randid).upnbrs.size()
 && all.nodes.get(randid).downnbrs.size() == all_random.nodes.get(randid).downnbrs.size()) {
 all_mirnas.remove(randid);
 }

 if (all_random_degree1 < all_degree1) {
 break;
 }
 }
 }

 if (it.int1.type.toLowerCase().equals("compound")) {
 System.out.println("----------4");
 while (true) {
 randint = random_number.nextInt(all_compounds.size());
 randid = all_compounds.get(randint);
 System.out.println("----------44 " + randint + " " + randid);
 if (dir.equals("indirect") || type.equals("kegg_reaction_reversible")) {
 all_degree1 = all.nodes.get(randid).revnbrs.size();
 all_random_degree1 = all_random.nodes.get(randid).revnbrs.size();
 } else {
 all_degree1 = all.nodes.get(randid).downnbrs.size();
 all_random_degree1 = all_random.nodes.get(randid).downnbrs.size();
 }
 System.out.println(all_degree1 + " " + all_random_degree1);
 if (all.nodes.get(randid).revnbrs.size() == all_random.nodes.get(randid).revnbrs.size()
 && all.nodes.get(randid).upnbrs.size() == all_random.nodes.get(randid).upnbrs.size()
 && all.nodes.get(randid).downnbrs.size() == all_random.nodes.get(randid).downnbrs.size()) {
 all_compounds.remove(randid);
 }
 if (all_random_degree1 < all_degree1) {
 break;
 }

 }
 }

 if (it.int1.type.toLowerCase().equals("glycan")) {
 System.out.println("----------5");
 while (true) {
 randint = random_number.nextInt(all_glycans.size());
 randid = all_glycans.get(randint);
 if (dir.equals("indirect") || type.equals("kegg_reaction_reversible")) {
 all_degree1 = all.nodes.get(randid).revnbrs.size();
 all_random_degree1 = all_random.nodes.get(randid).revnbrs.size();
 } else {
 all_degree1 = all.nodes.get(randid).downnbrs.size();
 all_random_degree1 = all_random.nodes.get(randid).downnbrs.size();
 }
 if (all.nodes.get(randid).revnbrs.size() == all_random.nodes.get(randid).revnbrs.size()
 && all.nodes.get(randid).upnbrs.size() == all_random.nodes.get(randid).upnbrs.size()
 && all.nodes.get(randid).downnbrs.size() == all_random.nodes.get(randid).downnbrs.size()) {
 all_glycans.remove(randid);
 }
 if (all_random_degree1 < all_degree1) {
 break;
 }
 }
 }

 if (it.int1.type.toLowerCase().equals("complex")) {
 System.out.println("----------6");
 while (true) {
 randint = random_number.nextInt(all_complexes.size());
 randid = all_complexes.get(randint);
 if (dir.equals("indirect") || type.equals("kegg_reaction_reversible")) {
 all_degree1 = all.nodes.get(randid).revnbrs.size();
 all_random_degree1 = all_random.nodes.get(randid).revnbrs.size();
 } else {
 all_degree1 = all.nodes.get(randid).downnbrs.size();
 all_random_degree1 = all_random.nodes.get(randid).downnbrs.size();
 }
 if (all.nodes.get(randid).revnbrs.size() == all_random.nodes.get(randid).revnbrs.size()
 && all.nodes.get(randid).upnbrs.size() == all_random.nodes.get(randid).upnbrs.size()
 && all.nodes.get(randid).downnbrs.size() == all_random.nodes.get(randid).downnbrs.size()) {
 all_complexes.remove(randid);
 }
 if (all_random_degree1 < all_degree1) {
 break;
 }
 }
 }

 if (it.int1.type.toLowerCase().equals("stimulus") || it.int1.type.toLowerCase().equals("chemical")
 || it.int1.type.toLowerCase().equals("smallmolecule")
 || it.int1.type.toLowerCase().equals("phenotype")
 || it.int1.type.toLowerCase().equals("proteinfamily")) {
 continue;
 }
 */
 /*int third_q = (int) Math.round(degree_distribution.length * 75.0 / 100.0);
 int first_q = (int) Math.round(degree_distribution.length * 25.0 / 100.0);
 double iqr = degree_distribution[third_q] - degree_distribution[first_q];
 double step = 2 * iqr * Math.pow(degree_distribution.length, -0.33);
 int number_of_bins = (int) Math.round((degree_distribution[degree_distribution.length - 1] - degree_distribution[0]) * (1.0 / step));
 System.out.println("\n" + degree_distribution[third_q] + " " + degree_distribution[first_q]);
 System.out.println(1 / step);
 System.out.println(step); */
 /*

 public void permute_phenotype_label(Network all, NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, String folder, String prefix, int count, List<String> hit_list, Map<String, String[]> hugo_by_id, String fplayers) throws IOException {

 System.out.println("+++++++++++++build_paths_for_random_hit_lists++++++++++++");
 java.util.Random randomGenerator = new java.util.Random();
 List<String> all_hugo_ids = new ArrayList();
 List<String> all_mirnas_ids = new ArrayList();
 int count_list_size = hit_list.size();
 int randomInt, degree;
 int max_degree = 0;
 int bin = 0;
 int index = 0;
 int number_of_proteins = 0;
 for (Node n : all.nodes.values()) {
 if (n.type.equals("protein")) {
 degree = n.downnbrs.size() + n.upnbrs.size() + n.revnbrs.size();
 if (degree > 0) {
 number_of_proteins++;
 }
 }
 }
 int[] degree_distribution = new int[number_of_proteins];
 for (Node n : all.nodes.values()) {
 if (n.type.equals("protein")) {
 degree = n.downnbrs.size() + n.upnbrs.size() + n.revnbrs.size();
 if (degree > 0) {
 degree_distribution[index] = degree;
 index++;
 if (max_degree < degree) {
 max_degree = degree;
 }
 }
 }
 }

 Arrays.sort(degree_distribution);
 int prev = degree_distribution[0];
 int number = 1;
 for (int i = 1; i < degree_distribution.length; i++) {
 if (prev != degree_distribution[i]) {
 System.out.print(prev + "-" + number + " ");
 prev = degree_distribution[i];
 number = 1;
 } else {
 number++;
 }
 }
 System.out.print(prev + "-" + number + "\n");

 int number_of_bins = 5;
 ArrayList<String>[] bins = new ArrayList[number_of_bins];
 for (bin = 0; bin < number_of_bins; bin++) {
 bins[bin] = new ArrayList();
 //System.out.println(bin + " " + bins[bin].size());
 }
 bin = 0;
 String id, suf;
 for (Node n : all.nodes.values()) {
 if (n.type.equals("protein")) {
 degree = n.downnbrs.size() + n.revnbrs.size() + n.upnbrs.size();
 if (degree > 0) {
 try {
 bin = get_bin(degree);
 if (bin == number_of_bins) {
 bin = number_of_bins;
 }
 bins[bin].add(hugo_by_id.get(n.id.substring(1))[1]);
 //all_hugo_ids.add(hugo_by_id.get(n.id.substring(1))[1]);
 } catch (NullPointerException e) {

 }
 }
 }
 if (n.type.toLowerCase().equals("miRNA".toLowerCase())) {
 id = n.id.replace("-3p", "").replace("-5p", "").replace("*", "");
 if (id.charAt(id.length() - 2) == '-' && Character.isDigit(id.charAt(id.length() - 1)) && Character.isDigit(id.charAt(id.length() - 3))) {
 suf = id.substring(id.length() - 2);
 id = id.substring(0, id.length() - 2);
 } else {
 suf = "";
 }
 if (id.contains("let-")) {
 all_mirnas_ids.add("MIR" + id.replace("-", "").substring(3).toUpperCase() + suf);
 System.out.println(id + " MIR" + id.replace("-", "").substring(3).toUpperCase());
 } else {
 all_mirnas_ids.add(id.replace("-", "").substring(3).toUpperCase() + suf);
 System.out.println(id + " " + id.replace("-", "").substring(3).toUpperCase());
 }
 /*
 System.out.println("\n\n\n" + n.id + " " +n.db_flag);
                
 for (String [] s: n.ids){
 for (int i=0; i< s.length; i++){
 System.out.println(s[i]);
 }
 }
 try {
 all_mirnas_ids.add(n.id);
 } catch (NullPointerException e) {
 } */ /*
 }
 }

 for (bin = 0; bin < number_of_bins; bin++) {
 System.out.println(bin + " " + bins[bin].size());
 }
 String next_gene;
 String flag;
 bin = 0;
 BufferedWriter random_hg_file = new BufferedWriter(new FileWriter(folder + "RandomHitGenes.txt"));
 BufferedWriter random_hg_file_lists = new BufferedWriter(new FileWriter(folder + "RandomHitGenes_lists.txt"));
 List<Integer>[] already_encountered_p = new ArrayList[10];
 List<Integer> already_encountered_m = new ArrayList();
 for (bin = 0; bin < number_of_bins; bin++) {
 already_encountered_p[bin] = new ArrayList();
 }
 for (int j = count; j >= 0; j--) {
 random_hg_file_lists.write("HitList" + j + "\n");
 for (int i = count_list_size - 1; i >= 0; i--) {

 if (!hit_list.get(i).startsWith("MIR")) {
 try {
 int h = all.nodes.get(hit_list.get(i)).downnbrs.size();
 } catch (NullPointerException e) {
 System.out.println(hit_list.get(i));
 }
 degree = all.nodes.get(hit_list.get(i)).downnbrs.size() + all.nodes.get(hit_list.get(i)).upnbrs.size() + all.nodes.get(hit_list.get(i)).revnbrs.size();
 bin = get_bin(degree);
 if (bin == number_of_bins) {
 bin = number_of_bins;
 }
 randomInt = randomGenerator.nextInt(bins[bin].size());
 next_gene = bins[bin].get(randomInt);
 flag = "p";
 } else {
 randomInt = randomGenerator.nextInt(all_mirnas_ids.size());
 next_gene = all_mirnas_ids.get(randomInt);
 flag = "m";
 }
 if ((flag.equals("p") && already_encountered_p[bin].contains(randomInt)) || (flag.equals("m") && already_encountered_m.contains(randomInt))) {
 random_hg_file_lists.write(next_gene + "\n");
 } else {
 random_hg_file_lists.write(next_gene + "\n");
 random_hg_file.write(next_gene + "\n");
 if (flag.equals("p")) {
 already_encountered_p[bin].add(randomInt);
 } else {
 already_encountered_m.add(randomInt);
 }
 }
 }
 }
 random_hg_file.close();
 random_hg_file_lists.close();

 nutils.load_hitlist_and_finalimpl(folder + "RandomHitGenes.txt", fplayers, dbutils.hugo);
 nutils.find_pathway_for_list_BF_algorithm(all, dbutils.hugo, 5, folder + "foundDF.txt", prefix);
 putils.find_the_shortest_paths(folder + "foundDF.txt", "_top_ranked", all, nutils.hg, nutils.fpl, 0, 0, 0, 0);

 }



 */
