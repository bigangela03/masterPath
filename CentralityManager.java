/*

 */
package masterPATH;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * CentralityManager class contains methods to calculate centrality scores
 *
 * @author Natalia Rubanova
 */
public class CentralityManager {

    /**
     * Find nodes and calculate centrality scores
     *
     * @param all Network
     * @param hugo_by_id Map with HGNC ids
     * @param nutils Link to NetworkManager object
     * @param f_shortest_pathways File with pathways
     * @param output_file Output file
     */
    public void calculate_centrality_scores_for_nodes(Network all, Map<String, String[]> hugo_by_id, NetworkManager nutils, String f_shortest_pathways, String output_file) throws FileNotFoundException, IOException {
        BufferedReader rd = new BufferedReader(new FileReader(f_shortest_pathways));
        BufferedWriter wr = new BufferedWriter(new FileWriter(output_file));
        String s;
        String id1;
        String id2;
        String[] ss;
        int i;
        List<String> tmp_hub_stat;
        Map<String, Integer> hub_stat = new HashMap();
        while ((s = rd.readLine()) != null) {
            tmp_hub_stat = new ArrayList();
            ss = s.split("\t");
            for (i = 3; i < ss.length - 2; i++) {
                id1 = all.interactions.get(ss[i]).int1.id;
                id2 = all.interactions.get(ss[i]).int2.id;
                if (!tmp_hub_stat.contains(id1)) {
                    tmp_hub_stat.add(id1);
                }
                if (!tmp_hub_stat.contains(id2)) {
                    tmp_hub_stat.add(id2);
                }
            }
            for (String t : tmp_hub_stat) {
                if (hub_stat.containsKey(t)) {
                    i = hub_stat.get(t);
                    hub_stat.put(t, i + 1);
                } else {
                    hub_stat.put(t, 1);
                }
            }
        }
        String s2;
        for (String t : hub_stat.keySet()) {
            if (t.startsWith("p")) {
                try {
                    s2 = "p" + hugo_by_id.get(t.substring(1))[1];
                } catch (NullPointerException e) {
                    s2 = t;
                }
            } else {
                try {
                    s2 = hugo_by_id.get(t)[1];
                } catch (NullPointerException e) {
                    s2 = t;
                }
            }
            wr.write(t + "\t" + s2 + "\t" + hub_stat.get(t) + "\n");
        }
        wr.close();
        rd.close();
    }
    
    /**
     * Find short linear paths in a mixed direct/indirect subnetwork and
     * calculate centrality scores
     *
     * @param all Network
     * @param inf List of pathways to create subnetwork
     * @param outf Output file
     * @param min_length minimum length of a path
     * @param max_length maximum length of a path
     */
    public void calculate_centrality_scores_for_paths(Network all, String inf, String outf, int min_length, int max_length) throws FileNotFoundException, IOException {
        System.out.println("+++++++++++++get_centrality_scores_for_paths++++++++++++");
        BufferedReader rd = new BufferedReader(new FileReader(inf));
        Map<S3, List<String>> all_paths_s3_id = new HashMap();
        Map<S3, List<Map<String, I1S1>>> common_parts_by_hg_fp_type = new HashMap();
        List<Map<String, I1S1>> tmp;
        Map<String, I2S2> common_parts = new HashMap();
        String s;
        String[] ss;
        S3 hg_fp_type_id;
        List<String> tmp_list_of_paths;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            try {
                hg_fp_type_id = new S3(ss[1], ss[ss.length - 1], returnPathType(all, s));
            } catch (ArrayIndexOutOfBoundsException e) {
                System.out.println(s);
                break;
            }
            if (!all_paths_s3_id.containsKey(hg_fp_type_id)) {
                tmp_list_of_paths = new ArrayList();
                tmp_list_of_paths.add(s);
                all_paths_s3_id.put(hg_fp_type_id, tmp_list_of_paths);
            } else {
                tmp_list_of_paths = all_paths_s3_id.get(hg_fp_type_id);
                tmp_list_of_paths.add(s);
                all_paths_s3_id.put(hg_fp_type_id, tmp_list_of_paths);
            }
        }

        for (S3 hg_fp_type : all_paths_s3_id.keySet()) {
            if (common_parts_by_hg_fp_type.containsKey(hg_fp_type)) {
                tmp = common_parts_by_hg_fp_type.get(hg_fp_type);
                tmp.add(find_unique_paths_for_list_of_pathways(max_length, min_length, all_paths_s3_id.get(hg_fp_type)));
                common_parts_by_hg_fp_type.put(hg_fp_type, tmp);
            } else {
                tmp = new ArrayList();
                tmp.add(find_unique_paths_for_list_of_pathways(max_length, min_length, all_paths_s3_id.get(hg_fp_type)));
                common_parts_by_hg_fp_type.put(hg_fp_type, tmp);
            }
        }

        for (S3 hg_fp_type : common_parts_by_hg_fp_type.keySet()) {
            for (Map<String, I1S1> part_type_pathid_map : common_parts_by_hg_fp_type.get(hg_fp_type)) {
                for (String t : part_type_pathid_map.keySet()) {
                    if (common_parts.containsKey(t)) {
                        common_parts.get(t).i1++;
                        //System.out.println(t);
                        if (common_parts.get(t).i2 != part_type_pathid_map.get(t).i1 && common_parts.get(t).i2 != 0) {
                            //common_parts.get(t).i2 = common_parts.get(t).i2 + 10;
                        }
                        if (common_parts.get(t).i2 != part_type_pathid_map.get(t).i1 && common_parts.get(t).i2 == 0) {
                            //common_parts.get(t).i2 = common_parts.get(t).i2 * 10;
                        }
                        common_parts.get(t).s1 = common_parts.get(t).s1 + ";" + hg_fp_type.s1 + "-" + hg_fp_type.s2 + "_" + hg_fp_type.s3;
                        common_parts.get(t).s2 = common_parts.get(t).s2 + ";" + part_type_pathid_map.get(t).s1;
                    } else {
                        common_parts.put(t, new I2S2(1, part_type_pathid_map.get(t).i1, hg_fp_type.s1 + "-" + hg_fp_type.s2 + "_" + hg_fp_type.s3, part_type_pathid_map.get(t).s1));
                    }
                }
            }
        }

        BufferedWriter out = new BufferedWriter(new FileWriter(outf));
        for (String t : common_parts.keySet()) {
            if (common_parts.get(t).i1 > 1) {
                ss = t.split("\t");
                out.write(t);
                for (int i = max_length; i >= ss.length; i--) {
                    out.write("\t");
                }
                out.write(common_parts.get(t).i2 + "\t" + common_parts.get(t).i1 + "\t" + common_parts.get(t).s1 + "\t" + common_parts.get(t).s2);
                out.newLine();
            }
        }
        out.close();
        rd.close();
    }

   

    /**
     * Find short linear paths in a indirect subnetwork and calculate centrality
     * scores
     *
     * @param all Network
     * @param inf List of pathways to create subnetwork
     * @param outf Output file
     * @param min_length Minimum path length
     * @param max_length Maximum path length
     */
    public void calculate_centrality_scores_for_paths_ppi(Network all, String inf, String outf, int min_length, int max_length) throws FileNotFoundException, IOException {
        System.out.println("+++++++++++++get_overrepr_unique_paths++++++++++++");
        BufferedReader rd = new BufferedReader(new FileReader(inf));
        Map<S3, List<String>> all_paths_s3_id = new HashMap();
        Map<S3, List<Map<String, I1S1>>> common_parts_by_hg_fp_type = new HashMap();
        List<Map<String, I1S1>> tmp;
        Map<String, I2S2> common_parts = new HashMap();
        String s;
        String[] ss;
        S3 hg_fp_type_id;
        List<String> tmp_list_of_paths;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            hg_fp_type_id = new S3(ss[1], ss[ss.length - 1], returnPathType(all, s));
            if (!all_paths_s3_id.containsKey(hg_fp_type_id)) {
                tmp_list_of_paths = new ArrayList();
                tmp_list_of_paths.add(s);
                all_paths_s3_id.put(hg_fp_type_id, tmp_list_of_paths);
            } else {
                tmp_list_of_paths = all_paths_s3_id.get(hg_fp_type_id);
                tmp_list_of_paths.add(s);
                all_paths_s3_id.put(hg_fp_type_id, tmp_list_of_paths);
            }
        }

        for (S3 hg_fp_type : all_paths_s3_id.keySet()) {
            if (common_parts_by_hg_fp_type.containsKey(hg_fp_type)) {
                tmp = common_parts_by_hg_fp_type.get(hg_fp_type);
                tmp.add(find_unique_paths_for_list_of_pathways_ppi(max_length, min_length, all_paths_s3_id.get(hg_fp_type)));
                common_parts_by_hg_fp_type.put(hg_fp_type, tmp);
            } else {
                tmp = new ArrayList();
                tmp.add(find_unique_paths_for_list_of_pathways_ppi(max_length, min_length, all_paths_s3_id.get(hg_fp_type)));
                common_parts_by_hg_fp_type.put(hg_fp_type, tmp);
            }
        }
        String symmetric;
        for (S3 hg_fp_type : common_parts_by_hg_fp_type.keySet()) {
            for (Map<String, I1S1> part_type_pathid_map : common_parts_by_hg_fp_type.get(hg_fp_type)) {
                for (String t : part_type_pathid_map.keySet()) {
                    symmetric = return_symmetric_path(t);
                    if (common_parts.containsKey(t)) {
                        common_parts.get(t).i1++;
                        //System.out.println(t);
                        if (common_parts.get(t).i2 != part_type_pathid_map.get(t).i1 && common_parts.get(t).i2 != 0) {
                            //common_parts.get(t).i2 = common_parts.get(t).i2 + 10;
                        }
                        if (common_parts.get(t).i2 != part_type_pathid_map.get(t).i1 && common_parts.get(t).i2 == 0) {
                            //common_parts.get(t).i2 = common_parts.get(t).i2 * 10;
                        }
                        common_parts.get(t).s1 = common_parts.get(t).s1 + ";" + hg_fp_type.s1 + "-" + hg_fp_type.s2 + "_" + hg_fp_type.s3;
                        common_parts.get(t).s2 = common_parts.get(t).s2 + ";" + part_type_pathid_map.get(t).s1;
                    } else if (common_parts.containsKey(symmetric)) {
                        common_parts.get(symmetric).i1++;
                        common_parts.get(symmetric).s1 = common_parts.get(symmetric).s1 + ";" + hg_fp_type.s1 + "-" + hg_fp_type.s2 + "_" + hg_fp_type.s3;
                        common_parts.get(symmetric).s2 = common_parts.get(symmetric).s2 + ";" + part_type_pathid_map.get(t).s1;
                    } else {
                        common_parts.put(t, new I2S2(1, part_type_pathid_map.get(t).i1, hg_fp_type.s1 + "-" + hg_fp_type.s2 + "_" + hg_fp_type.s3, part_type_pathid_map.get(t).s1));
                    }
                }
            }
        }

        BufferedWriter out = new BufferedWriter(new FileWriter(outf));
        for (String t : common_parts.keySet()) {
            ss = t.split("\t");
            out.write(t);
            for (int i = max_length; i >= ss.length; i--) {
                out.write("\t");
            }
            out.write(common_parts.get(t).i2 + "\t" + common_parts.get(t).i1 + "\t" + common_parts.get(t).s1 + "\t" + common_parts.get(t).s2);
            out.newLine();
        }
        out.close();
        rd.close();
    }
    /**
     * Add approved gene symbols to file with paths and centrality scores
     *
     * @param inf Input file with paths
     * @param outf Output file
     * @param all Network all
     * @param hugo_by_id Map of approved gene symbols and gene IDs from
     * DBManager class
     */
    public void add_hgnc_symbols_to_paths(String inf, String outf, Network all, Map<String, String[]> hugo_by_id) throws IOException {
        System.out.println("+++++++++++++put_names_for_overreprs_unique++++++++++++");
        BufferedReader rd = new BufferedReader(new FileReader(inf));
        BufferedWriter wr = new BufferedWriter(new FileWriter(outf));
        String s, n1, n2, scon1 = "", scon2 = "", scon, sn1, sn2;
        String[] ss;
        n1 = "";
        n2 = "";
        int i = 0;

        String paths, path_names = "", t;
        String[] pths, tmp1, tmp2;
        int len = 0, con = 0, con1 = 0, con2 = 0;
        while ((s = rd.readLine()) != null) {
            con = 0;
            scon = "";
            ss = s.split("\t");
            path_names = "";
            i = 0;
            //for (int i = 0; i < ss.length - 2; i++) {
            paths = ss[ss.length - 2];
            pths = paths.split(";");
            String s1, s2;
            for (String p : pths) {
                tmp1 = p.split("_");
                tmp2 = tmp1[0].split("-");
                try {
                    if (tmp2[0].startsWith("p")) {
                        s1 = "p" + hugo_by_id.get(tmp2[0].substring(1))[1];
                    } else {
                        s1 = hugo_by_id.get(tmp2[0])[1];
                    }
                    if (tmp2[1].startsWith("p")) {
                        s2 = "p" + hugo_by_id.get(tmp2[1].substring(1))[1];
                    } else {
                        s2 = hugo_by_id.get(tmp2[1])[1];
                    }
                    path_names = path_names + ";" + s1 + "-" + s2 + "_" + tmp1[1];
                } catch (NullPointerException e) {
                    System.out.println(p);
                    break;
                }
            }
            t = "";
            len = 0;
            while ((i < ss.length - 4) && !ss[i].equals("")) {
                len++;
                try {
                    if (all.interactions.get(ss[i]).int1.id.startsWith("hsa-")) {
                        n1 = all.interactions.get(ss[i]).int1.id;
                    } else if (all.interactions.get(ss[i]).int1.id.startsWith("p")) {
                        n1 = "p" + hugo_by_id.get(all.interactions.get(ss[i]).int1.id.substring(1))[1];
                    } else {
                        n1 = hugo_by_id.get(all.interactions.get(ss[i]).int1.id)[1];
                    }
                } catch (NullPointerException e) {

                    try {
                        n1 = all.interactions.get(ss[i]).int1.id;
                    } catch (NullPointerException ee) {
                        n1 = "-";
                        System.out.println("unique names exception " + ss[i]);
                    }
                }
                try {
                    if (all.interactions.get(ss[i]).int2.id.startsWith("hsa-")) {
                        n2 = all.interactions.get(ss[i]).int2.id;
                    } else if (all.interactions.get(ss[i]).int2.id.startsWith("p")) {
                        n2 = "p" + hugo_by_id.get(all.interactions.get(ss[i]).int2.id.substring(1))[1];
                    } else {
                        n2 = hugo_by_id.get(all.interactions.get(ss[i]).int2.id)[1];
                    }
                } catch (NullPointerException e) {
                    // System.out.println(all.interactions.get(ss[i]).int2.id);
                    n2 = all.interactions.get(ss[i]).int2.id;
                    //n2 = "-";//all.interactions.get(ss[i]).int1.id;
                    //System.out.println("B\n" + ss[i]);
                    //System.out.println(s);
                }

                // wr.write("\t" + n1 + "-" + n2);
                t = t + "\t" + n1 + "-" + n2;
                i++;
            }

            for (int j = 0; j < len; j++) {
                try {
                    con1 = (all.interactions.get(ss[j]).int1.downnbrs.size() + all.interactions.get(ss[j]).int1.revnbrs.size()) * (all.interactions.get(ss[j]).int1.upnbrs.size() + all.interactions.get(ss[j]).int1.revnbrs.size());
                    con2 = (all.interactions.get(ss[j]).int2.downnbrs.size() + all.interactions.get(ss[j]).int2.revnbrs.size()) * (all.interactions.get(ss[j]).int2.upnbrs.size() + all.interactions.get(ss[j]).int2.revnbrs.size());
                    scon1 = "" + (all.interactions.get(ss[j]).int1.downnbrs.size() + all.interactions.get(ss[j]).int1.revnbrs.size()) + "*" + (all.interactions.get(ss[j]).int1.upnbrs.size() + all.interactions.get(ss[j]).int1.revnbrs.size()) + " " + all.interactions.get(ss[j]).int1.id + " " + returnName(all, ss[j], hugo_by_id, 1);
                    scon2 = "" + (all.interactions.get(ss[j]).int2.downnbrs.size() + all.interactions.get(ss[j]).int2.revnbrs.size()) + "*" + (all.interactions.get(ss[j]).int2.upnbrs.size() + all.interactions.get(ss[j]).int2.revnbrs.size()) + " " + all.interactions.get(ss[j]).int2.id + " " + returnName(all, ss[j], hugo_by_id, 2);
                    //System.out.println(scon1 + "\n" +scon2);
                } catch (NullPointerException e) {
                    System.out.println(j + "---" + ss[j]);
                }
                if (con1 > con) {
                    con = con1;
                    scon = scon1;
                }
                if (con2 > con) {
                    con = con2;
                    scon = scon2;
                }
                //   System.out.println(scon);
            }

            wr.write(len + "\t" + s + "\t" + con + "\t" + scon + "\t" + path_names + t);
            // System.out.println(scon);
            wr.newLine();
        }
        rd.close();
        wr.close();
    }
    /**
     * Add scores and p-values from GW screen to file with nodes and centrality scores
     *
     * @param nodes_file File with nodes and centrality scores
     * @param gw_file File with screening data
     * @param conv_table Conversion table for genes ids
     * @param out Output file
     */
    public void add_experimental_scores_and_pvalues_to_nodes(String nodes_file, String gw_file, String conv_table, String out) throws FileNotFoundException, IOException {
        BufferedReader rd_cp = new BufferedReader(new FileReader(nodes_file));
        BufferedWriter wr = new BufferedWriter(new FileWriter(out));
        String s;
        String[] ss;
        Map<String, Float> p_values = new HashMap();
        Map<String, Float> fold_changes = new HashMap();
        Map<String, Float[]> table ;

        table = load_screening_data(gw_file, conv_table);

        for (String q : table.keySet()) {
            p_values.put(q, table.get(q)[1]);
            fold_changes.put(q, table.get(q)[0]);
        }
        String name;
        while ((s = rd_cp.readLine()) != null) {
            ss = s.split("\t");
            if (Float.parseFloat(ss[3])>=0.05){
                continue;
            }
            if (ss[1].startsWith("p")) {
                name = ss[1].substring(1);
            } else if (ss[1].startsWith("hsa-")) {
                name = ss[1];
            } else {
                name = ss[1];
            }

            if (p_values.containsKey(name)) {
                wr.write(s + "\t" + fold_changes.get(name) + "\t" + p_values.get(name)+"\n");
            } else {
                wr.write(s + "\t-\t-\n");
            }

        }
        rd_cp.close();
        wr.close();
    }

    /**
     * Add aggregated p-values from GW screen to file with paths and centrality
     * scores
     * @param nodes_file File with nodes and centrality scores
     * @param gw_file File with screening data
     * @param conv_table Conversion table for genes ids
     * @param out Output file
     */
    public void add_experimnetal_scores_and_aggregated_pvalues_to_paths(String unique_over, String gw_file, String conv_table, String out) throws FileNotFoundException, IOException {
        BufferedReader rd_cp = new BufferedReader(new FileReader(unique_over));
        //BufferedReader rd_pv = new BufferedReader(new FileReader(gw_file));
        BufferedWriter wr = new BufferedWriter(new FileWriter(out));
        int pos = 12;
        String s;
        String[] ss, tt;
        int len, number;
        double pvalue, tfct, tfct_test, foldchange;
        double[] tfct_permuted = new double[1000];
        Map<String, List<Float>> permutation_table = new HashMap();
        Map<String, Float> p_values = new HashMap();
        Map<String, Float> fold_changes = new HashMap();
        List<String> genes, genes_new;
        permutation_table = return_permutation_distribution(gw_file, conv_table, 10000);
        Map<String, Float[]> table = new HashMap();

        table = load_screening_data(gw_file, conv_table);

        for (String q : table.keySet()) {
            p_values.put(q, table.get(q)[1]);
            fold_changes.put(q, table.get(q)[0]);
        }

        while ((s = rd_cp.readLine()) != null) {

            genes = new ArrayList();
            genes_new = new ArrayList();
            ss = s.split("\t");

            len = Integer.valueOf(ss[0]);
            if (Float.parseFloat(ss[ss.length - 1]) >= 0.05 || Integer.parseInt(ss[6]) <= 3 || Integer.parseInt(ss[0]) < 2) {
                continue;
            }
            wr.write(s + "\t");
            for (int i = 0; i < len; i++) {
                if (!ss[pos + i].contains("hsa-")) {
                    tt = ss[pos + i].split("-");
                    if (tt[0].startsWith("p")) {
                        tt[0] = tt[0].substring(1);
                    }
                    if (tt[1].startsWith("p")) {
                        tt[1] = tt[1].substring(1);
                    }
                    if (!genes.contains(tt[0])) {
                        genes.add(tt[0]);
                    }
                    if (!genes.contains(tt[1])) {
                        genes.add(tt[1]);
                    }
                } else {
                    tt = ss[pos + i].split("-");
                    //System.out.print (ss[pos+i]  );
                    for (int j = 0; j < tt.length; j++) {
                        if (ss[pos + i].startsWith("hsa-") && Character.isUpperCase(tt[j].charAt(0)) && !genes.contains(tt[j])) {
                            genes.add(tt[j]);
                        }
                        try {
                            if (ss[pos + i].startsWith("p") && Character.isUpperCase(tt[0].charAt(1)) && !genes.contains(tt[j])) {
                                genes.add(tt[j].substring(1));
                            }
                        } catch (StringIndexOutOfBoundsException e) {
                            System.out.println(tt[j] + " " + pos + " " + i + " " + s);
                            break;
                        }
                    }
                }
                // System.out.println(tt[0] + " " + tt[1]);
            }
            tfct = 0.0;
            foldchange = 0.0;
            for (String gene : genes) {
                if (p_values.containsKey(gene)) {
                    genes_new.add(gene);
                    //System.out.println(p_values.get(gene));
                    tfct = tfct + (Math.log(p_values.get(gene)) * 1.0) / Math.log(2);
                    foldchange = foldchange + fold_changes.get(gene);
                } else {
                    //System.out.println(gene + " ");
                }
            }
            tfct = -2 * tfct;
            foldchange = foldchange / genes.size();
            number = 0;
            if (genes.size() == 1) {
                // System.out.println(s);
            }
            // System.out.println(genes.size() + " " + genes_new.size() +" " + tfct);
            for (int i = 0; i < 1000; i++) {
                tfct_test = 0.0;
                for (String gene : genes_new) {
                    tfct_test = tfct_test + (Math.log(permutation_table.get(gene).get(i)) * 1.0) / Math.log(2);
                }
                tfct_test = -2 * tfct_test;
                if (tfct_test >= tfct) {
                    number++;
                }
                //System.out.print(" " + tfct_test);
                tfct_permuted[i] = tfct_test;
            }
            //System.out.println();
            pvalue = (number * 1.0) / 1000;

            wr.write(foldchange + "\t" + pvalue + "\n");

        }
        wr.close();
        rd_cp.close();
        // rd_pv.close();
    }

    /**
     * Add information about degree to file with paths and centrality scores
     *
     * @param unique_over File with paths
     * @param path_con File with pathway connectivity
     * @param out Output file
     */
    public void calculate_paths_degree(String unique_over, String path_con, String out) throws FileNotFoundException, IOException {
        BufferedReader rd_uo = new BufferedReader(new FileReader(unique_over));
        BufferedReader rd_pc = new BufferedReader(new FileReader(path_con));
        BufferedWriter wr = new BufferedWriter(new FileWriter(out));
        String s;
        String[] ss, tt;
        int len;
        Map<String, String> pthw_con = new HashMap();
        while ((s = rd_pc.readLine()) != null) {
            ss = s.split("\t");
            pthw_con.put(ss[0], ss[ss.length - 1]);
        }
        while ((s = rd_uo.readLine()) != null) {
            ss = s.split("\t");
            len = Integer.valueOf(ss[0]);
            for (int i = 0; i <= len + 4; i++) {
                // wr.write(ss[i] + "\t");
            }
            wr.write(s + "\t");
            tt = ss[8].split(";");
            // System.out.println(ss[8]);
            for (int i = 0; i < tt.length - 1; i++) {
                wr.write(pthw_con.get(tt[i]) + ";");
                //System.out.println(tt[i] + pthw_con.get(tt[i]));
            }
            wr.write(pthw_con.get(tt[tt.length - 1]) + "\n");
            for (int i = len + 5; i < ss.length - 2; i++) {
                //wr.write(ss[i] + "\t");
            }
            //wr.write(ss[ss.length-1]+"\n");
        }
        wr.close();
        rd_pc.close();
        rd_uo.close();
    }
    
    /**
     * Find all unique paths in a indirect subnetwork
     *
     * @param max_length Maximum length
     * @param min_length Minimum length
     * @param all_paths Pathways to create subnetwork
     */
    public Map<String, I1S1> find_unique_paths_for_list_of_pathways_ppi(int max_length, int min_length, List<String> all_paths) throws IOException {
        Map<String, I1S1> unique_parts = new HashMap();
        Map<String, I1S1> tmp_parts;
        List<String> paths = new ArrayList();
        String part, path_id;
        String[] ss;
        int min, max, max_end;
        paths.addAll(all_paths);
        for (String path : all_paths) {
            //System.out.println(path);
            ss = path.split("\t");
            if (ss.length - 3 < max_length) {
                max = ss.length - 3;
            } else {
                max = max_length;
            }
            min = min_length;
            paths.remove(path);
            path_id = ss[0];
            for (int i = 2; i <= ss.length - 1 - max; i++) {
                for (int j = min; j <= max; j++) {
                    part = ss[i];
                    for (int k = i + 1; k < i + j; k++) {
                        part = part + "\t" + ss[k];
                    }
                    tmp_parts = return_unique_path_ppi(paths, unique_parts, part, i, path_id);
                    unique_parts.putAll(tmp_parts);
                }
            }

            max_end = max;
            for (int i = ss.length - max; i <= ss.length - 1 - min; i++) { // <= -2????
                max_end = max_end - 1;
                for (int j = min; j <= max_end; j++) {
                    part = ss[i];
                    for (int k = i + 1; k < i + j; k++) {
                        part = part + "\t" + ss[k];
                    }
                    tmp_parts = return_unique_path_ppi(paths, unique_parts, part, i, path_id);
                    unique_parts.putAll(tmp_parts);
                }
            }
        }
        return unique_parts;
    }
     /**
     * Filter paths
     *
     * @param foundf File with paths
     * @param outname Output file
     * @param all Network
     * @param hg Hit genes
     * @param fpl Final Players
     * @param min_occ Minimum occurrence filter
     * @param mask Mask to include only paths with miRNAs
     */
    public void filter_paths(String foundf, String outname, Network all, Map<String, String> hg, Map<String, String> fpl, int min_occ, String mask) throws FileNotFoundException, IOException {
        System.out.println("+++++++++++++filter overrepresented paths++++++++++++");
        BufferedReader rd = new BufferedReader(new FileReader(foundf));
        BufferedWriter wr = new BufferedWriter(new FileWriter(foundf + outname));
        String s, hgene, fplr, mirna;
        String[] ss, tmp;
        boolean notmiddle;
        while ((s = rd.readLine()) != null) {
            notmiddle = false;
            ss = s.split("\t");
            if (mask.equals("occ")) {
                if (Integer.parseInt(ss[0]) >= min_occ) { //
                    wr.write(s);
                    wr.newLine();
                }
            } else if (mask.equals("mirna_midlle")) {
                if (!s.contains("hsa-")) {
                    continue;
                }
                hgene = ss[3];
                fplr = ss[4];
                for (int i = 4 + Integer.parseInt(ss[1]) + 1; i < ss.length; i++) {
                    if (!ss[i].contains("hsa-")) {
                        continue;
                    }
                    tmp = ss[i].split("-");
                    if (ss[i].startsWith("hsa-")) {
                        mirna = tmp[0] + "-" + tmp[1] + "-" + tmp[2];
                    } else {
                        mirna = tmp[1] + "-" + tmp[2] + "-" + tmp[3];
                    }
                    if (!mirna.equals(hgene) && !mirna.equals(fplr)) {
                        notmiddle = true;
                        break;
                    }
                }
                if (notmiddle) {
                    wr.write(s);
                    wr.newLine();
                }
            }
        }
        rd.close();
        wr.close();
    }
     /**
     * Load experimental data
     *
     * @param p_file Experimental data
     * @param n_file File for symbol conversion
     */
    public HashMap<String, Float[]> load_screening_data(String p_file, String n_file) throws IOException {
        BufferedReader rd = new BufferedReader(new FileReader(p_file));
        BufferedReader ct = new BufferedReader(new FileReader(n_file));

        String s;
        String[] ss;

        Float[] data = new Float[2];
        HashMap<String, Float[]> table = new HashMap();
        Map<String, String> conv_table = new HashMap();
        int index;

        while ((s = ct.readLine()) != null) {
            ss = s.split("\t");
            if (ss[1].equals("Unmatched")) {
                conv_table.put(ss[0], s + ss[0]);
                continue;
            }
            if (!conv_table.containsKey(ss[0])) {
                conv_table.put(ss[0], s);
            } else {
                if (conv_table.get(ss[0]).contains("Approved symbol")) {

                } else if (ss[1].equals("Approved symbol")) {
                    conv_table.put(ss[0], s);
                } else {
                    System.out.println(s);
                    System.out.println(conv_table.get(ss[0]));
                }
            }
        }

        String apname = "";
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            if (ss[3].contains("//")) {
                ss[3] = ss[3].split("//")[0];
            }
            if (!ss[3].equals("---")) {
                data = new Float[2];
                try {
                    apname = conv_table.get(ss[3]).split("\t")[2];
                } catch (NullPointerException e) {
                    System.out.println("----" + ss[3]);
                }
                if (!table.containsKey(apname)) {
                    data[0] = Math.abs(Float.parseFloat(ss[0]));
                    data[1] = Float.parseFloat(ss[1]);
                    table.put(apname, data);
                } else {
                    if (Float.parseFloat(ss[1]) < table.get(apname)[1]) {
                        data[0] = Math.abs(Float.parseFloat(ss[0]));
                        data[1] = Float.parseFloat(ss[1]);
                        table.put(apname, data);
                        // table.put(apname, Float.parseFloat(ss[1]));
                    }
                }
            }
        }

        rd.close();
        ct.close();
        return table;
    }
     /**
     * Calculate empirical distribution for Fisher's statistic
     *
     * @param p_file Experimental data
     * @param n_file File for symbol conversion
     * @param number_of_permutations  Number of permutations
     */
    public HashMap<String, List<Float>> return_permutation_distribution(String p_file, String n_file, int number_of_permutations) throws FileNotFoundException, IOException {
        BufferedReader rd = new BufferedReader(new FileReader(p_file));
        BufferedReader ct = new BufferedReader(new FileReader(n_file));
        HashMap<String, List<Float>> permutation_table = new HashMap();
        String s;
        String[] ss;
        List<String> gnames = new ArrayList();
        List<Float> pvalues = new ArrayList();
        List<Float> test_list;
        HashMap<String, Float[]> table = new HashMap();

        int index;

        table = load_screening_data(p_file, n_file);

        for (String gene : table.keySet()) {
            test_list = new ArrayList();
            permutation_table.put(gene, test_list);
            pvalues.add(table.get(gene)[1]);
            gnames.add(gene);
        }
        String gene;
        for (int i = 0; i < number_of_permutations; i++) {
            Collections.shuffle(pvalues);
            for (index = 0; index < gnames.size(); index++) {
                gene = gnames.get(index);
                permutation_table.get(gene).add(pvalues.get(index));
            }
        }
        rd.close();
        ct.close();
        return permutation_table;
    }

    
    /**
     * Return symmetric pathway
     *
     * @param part Pathway
     */
    public static String return_symmetric_path(String part) {
        String[] ss;
        String[] tt;
        ss = part.split("\t");
        String symmetric = ss[ss.length - 1];
        for (int i = ss.length - 2; i >= 0; i--) {
            symmetric = symmetric + "\t" + ss[i];
        }
        return symmetric;
    }
     /**
     * Find all unique paths in a mixed direct/indirect subnetwork
     *
     * @param max_length Minimum length
     * @param min_length Maximum length
     * @param all_paths List of pathways to create subnetwork
     */
    public Map<String, I1S1> find_unique_paths_for_list_of_pathways(int max_length, int min_length, List<String> all_paths) throws IOException {
        Map<String, I1S1> unique_parts = new HashMap();
        Map<String, I1S1> tmp_parts;
        List<String> paths = new ArrayList();
        String part, path_id;
        String[] ss;
        int min, max, max_end;
        paths.addAll(all_paths);
        for (String path : all_paths) {
            //System.out.println(path);
            ss = path.split("\t");
            if (ss.length - 3 < max_length) {
                max = ss.length - 3;
            } else {
                max = max_length;
            }
            min = min_length;
            paths.remove(path);
            path_id = ss[0];
            for (int i = 2; i <= ss.length - 1 - max; i++) {
                for (int j = min; j <= max; j++) {
                    part = ss[i];
                    for (int k = i + 1; k < i + j; k++) {
                        part = part + "\t" + ss[k];
                    }
                    tmp_parts = return_unique_path(paths, unique_parts, part, i, path_id);
                    unique_parts.putAll(tmp_parts);
                }
            }

            max_end = max;
            for (int i = ss.length - max; i <= ss.length - 1 - min; i++) { // <= -2????
                max_end = max_end - 1;
                for (int j = min; j <= max_end; j++) {
                    part = ss[i];
                    for (int k = i + 1; k < i + j; k++) {
                        part = part + "\t" + ss[k];
                    }
                    tmp_parts = return_unique_path(paths, unique_parts, part, i, path_id);
                    unique_parts.putAll(tmp_parts);
                }
            }
        }
        return unique_parts;
    }

    /**
     * Compare two lists of nodes
     *
     * @param f1 File 1
     * @param f2 File 2
     */
    public void compare_two_lists_of_nodes(String f1, String f2) throws FileNotFoundException, IOException {
        BufferedReader rd1 = new BufferedReader(new FileReader("F:\\Dropbox\\_\u0420\u0430\u0431\u043e\u0442\u0430\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\output\\hubs_CC3_good_pvalue.txt"));
        BufferedReader rd2 = new BufferedReader(new FileReader("F:\\Dropbox\\_\u0420\u0430\u0431\u043e\u0442\u0430\\Programms\\DATA\\PhD\\Systems\\Case2_Control\\input\\HitGenes.txt"));
        String s;
        List<String> genes = new ArrayList();
        while ((s = rd1.readLine()) != null) {
            genes.add(s);
        }
        System.out.println("hubs " + genes.size());
        int count_hits = 0;
        int count_common = 0;
        while ((s = rd2.readLine()) != null) {
            count_hits++;
            if (genes.contains(s)) {
                System.out.println(s);
                count_common++;
            }
        }
        System.out.println("hits " + count_hits);
        System.out.println("common " + count_common);
    }
    
    /**
     * Check type of the pathway
     *
     * @param all Network
     * @param s pathway
     */
    public String returnPathType(Network all, String s) {
        String[] ss;
        String t;
        ss = s.split("\t");
        if (s.contains("KEGG")) {
            return "rn:";
        } else if (s.contains("MIRT0") || s.contains("transmir:") || s.contains("_mr")) {
            return "mirrna";
        } else if (s.contains("man")) {
            return "pg";
        } else {
            for (int i = 2; i < ss.length - 1; i++) {
                t = ss[i];
                // System.out.println(s+"\t"+t);

                try {
                    if (all.interactions.get(t).int1.type.equals("compound")
                            || all.interactions.get(t).int1.type.equals("glycan")
                            || all.interactions.get(t).int2.type.equals("compound")
                            || all.interactions.get(t).int1.type.equals("glycan")) {
                        return "rn:";
                    }
                    if (all.interactions.get(t).int1.type.equals("mirrna") || all.interactions.get(t).int2.type.equals("mirrna")) {
                        return "mirrna";
                    }
                    if (all.interactions.get(t).int1.type.equals("gene") || all.interactions.get(t).int2.type.equals("gene")) {
                        return "pg";
                    }
                } catch (NullPointerException e) {
                    System.out.println("----" + t + " " + s);
                    return "pg";
                }
            }
            return "pp";
        }

    }
    
    private Map<String, I1S1> return_unique_path(List<String> all_paths, Map<String, I1S1> unique_parts, String part, int offset_in_path, String path_id) {
        Map<String, I1S1> parts = new HashMap();
        I1S1 i1s1;
        if (!unique_parts.containsKey(part)) {
            int type = get_occurance_type(all_paths, part, offset_in_path);
            parts.put(part, new I1S1(type, path_id));
        } else {
            i1s1 = unique_parts.get(part);
            i1s1.s1 = i1s1.s1 + ";" + path_id; // automatically replacing in uniqe_parts because they have the same key
            parts.put(part, i1s1);
        }
        return parts;
    }

    private Map<String, I1S1> return_unique_path_ppi(List<String> all_paths, Map<String, I1S1> unique_parts, String part, int offset_in_path, String path_id) {
        Map<String, I1S1> parts = new HashMap();
        String symmetric = return_symmetric_path(part);
        I1S1 i1s1;
        if (!unique_parts.containsKey(part) && !unique_parts.containsKey(symmetric)) {
            int type = get_occurance_type(all_paths, part, offset_in_path);
            parts.put(part, new I1S1(type, path_id));
        } else if (unique_parts.containsKey(part)) {
            i1s1 = unique_parts.get(part);
            i1s1.s1 = i1s1.s1 + ";" + path_id; // automatically replacing in uniqe_parts because they have the same key
            parts.put(part, i1s1);
        } else if (unique_parts.containsKey(symmetric)) {
            i1s1 = unique_parts.get(symmetric);
            i1s1.s1 = i1s1.s1 + ";" + path_id; // automatically replacing in uniqe_parts because they have the same key
            parts.put(symmetric, i1s1);
        }
        return parts;
    }

    private int get_occurance_type(List<String> all_paths, String part, int offset_in_path) {
        String[] path_ss, part_ss;
        String tmp_part;
        int result = 0;
        for (String path : all_paths) {
            path_ss = path.split("\t");
            part_ss = part.split("\t");
            if (path_ss.length >= offset_in_path + part_ss.length) {
                tmp_part = path_ss[offset_in_path];
                for (int j = offset_in_path + 1; j < offset_in_path + part_ss.length; j++) {
                    tmp_part = tmp_part + "\t" + path_ss[j];
                }
                if (tmp_part.equals(part)) {
                    result = 1;
                }
            }
        }
        return result;
    }

   

    private String returnName(Network all, String interaction, Map<String, String[]> hugo_by_id, int i) {
        String n1;
        if (i == 1) {
            try {
                if (all.interactions.get(interaction).int1.id.startsWith("hsa-")) {
                    n1 = all.interactions.get(interaction).int1.id;
                } else if (all.interactions.get(interaction).int1.id.startsWith("p")) {
                    n1 = "p" + hugo_by_id.get(all.interactions.get(interaction).int1.id.substring(1))[1];
                } else {
                    n1 = hugo_by_id.get(all.interactions.get(interaction).int1.id)[1];
                }
            } catch (NullPointerException e) {

                try {
                    n1 = all.interactions.get(interaction).int1.id;
                } catch (NullPointerException ee) {
                    n1 = "-";
                    System.out.println("unique names exception " + interaction);
                }
            }
        } else {
            try {
                if (all.interactions.get(interaction).int2.id.startsWith("hsa-")) {
                    n1 = all.interactions.get(interaction).int2.id;
                } else if (all.interactions.get(interaction).int2.id.startsWith("p")) {
                    n1 = "p" + hugo_by_id.get(all.interactions.get(interaction).int2.id.substring(1))[1];
                } else {
                    n1 = hugo_by_id.get(all.interactions.get(interaction).int2.id)[1];
                }
            } catch (NullPointerException e) {

                try {
                    n1 = all.interactions.get(interaction).int2.id;
                } catch (NullPointerException ee) {
                    n1 = "-";
                    System.out.println("unique names exception " + interaction);
                }
            }
        }
        return n1;
    }


    

    private void find_overreprs(int max_length, int min_length, String inf, String outf) throws IOException {
        System.out.println("+++++++++++++find_overreprs++++++++++++");
        BufferedReader rd = new BufferedReader(new FileReader(inf));
        Map<Integer, List<Link>> links = new HashMap();
        Map<String, String[]> paths = new HashMap();
        String s;
        String[] ss, pth;
        //String tmp = "";
        int min, max, max_end;//, tmpi;
        int u = 0, count = 0;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            paths.put(ss[0], ss);
            //ii++;
        }
        for (String id : paths.keySet()) {
            count++;
            ss = paths.get(id);
            if (ss.length - 3 < max_length) {
                max = ss.length - 3;
            } else {
                max = max_length;
            }
            min = min_length;
            for (int offset = 2; offset <= ss.length - 1 - max; offset++) {
                for (int lnth = min; lnth <= max; lnth++) {
                    pth = new String[lnth];
                    u = 0;
                    for (int k = offset; k < offset + lnth; k++) {
                        //tmpa[u] = ss[k];
                        //pth[k - offset] = ss[k];
                        pth[u] = ss[k];
                        u++;
                    }
                    checkLinks(paths, links, pth, id, offset, lnth);
                }
            }
            max_end = max;
            for (int offset = ss.length - max; offset <= ss.length - 1 - min; offset++) { // <= -2????
                max_end = max_end - 1;
                for (int lnth = min; lnth <= max_end; lnth++) {
                    pth = new String[lnth];
                    u = 0;
                    for (int k = offset; k < offset + lnth; k++) {
                        //tmpa[u] = ss[k];
                        //pth[k - offset] = ss[k];
                        pth[u] = ss[k];
                        u++;
                    }
                    checkLinks(paths, links, pth, id, offset, lnth);
                }
            }
        }
        writeLinks(paths, links, outf);
        //System.out.println("links size " + links.size());
    }

    private void put_gene_symbols_for_paths(String inf, String outf, Network all, Map<String, String[]> hugo_by_id) throws IOException {
        System.out.println("+++++++++++++put_names_for_overreprs++++++++++++");
        BufferedReader rd = new BufferedReader(new FileReader(inf));
        BufferedWriter wr = new BufferedWriter(new FileWriter(outf));
        String s, n1, n2;
        String[] ss;
        n1 = "";
        n2 = "";
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            wr.write(s);
            for (int i = 5; i < ss.length; i++) {
                try {
                    if (all.interactions.get(ss[i]).int1.id.startsWith("hsa-")) {
                        n1 = all.interactions.get(ss[i]).int1.id;
                    } else if (all.interactions.get(ss[i]).int1.id.startsWith("p")) {
                        n1 = "p" + hugo_by_id.get(all.interactions.get(ss[i]).int1.id.substring(1))[1];
                    } else {
                        n1 = hugo_by_id.get(all.interactions.get(ss[i]).int1.id)[1];
                    }
                } catch (NullPointerException e) {
                    //System.out.println(all.interactions.get(ss[i]).int1.id);
                    try {
                        n1 = all.interactions.get(ss[i]).int1.id;

                    } catch (NullPointerException ee) {
                        System.out.println("A\n" + ss[i]);
                        System.out.println(s);
                        break;
                        //  n1 = "-";//all.interactions.get(ss[i]).int1.id;
                        //  System.out.println("A\n" + ss[i]);
                        //  System.out.println(s);
                    }
                }
                try {
                    if (all.interactions.get(ss[i]).int2.id.startsWith("hsa-")) {
                        n2 = all.interactions.get(ss[i]).int2.id;
                    } else if (all.interactions.get(ss[i]).int2.id.startsWith("p")) {
                        n2 = "p" + hugo_by_id.get(all.interactions.get(ss[i]).int2.id.substring(1))[1];
                    } else {
                        n2 = hugo_by_id.get(all.interactions.get(ss[i]).int2.id)[1];
                    }
                } catch (NullPointerException e) {
                    // System.out.println(all.interactions.get(ss[i]).int2.id);
                    n2 = all.interactions.get(ss[i]).int2.id;
                }
                wr.write("\t" + n1 + "-" + n2);
            }
            wr.newLine();
        }
        rd.close();
        wr.close();
    }

    private void find_paths_of_lenght1(String inf, String outf) throws FileNotFoundException, IOException {
        BufferedReader rd = new BufferedReader(new FileReader(inf));
        BufferedWriter wr = new BufferedWriter(new FileWriter(outf));
        Map<String, String[]> singles = new HashMap();
        String[] ss;
        String s;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            for (int i = 5; i < 5 + Integer.parseInt(ss[1]); i++) {
                if (singles.containsKey(ss[i])) {
                    if (Integer.parseInt(singles.get(ss[i])[0]) < Integer.parseInt(ss[0])) {
                        singles.put(ss[i], ss);
                    }
                } else {
                    singles.put(ss[i], ss);
                }
            }
        }
        for (String t : singles.keySet()) {
            wr.write(t + "\t" + singles.get(t)[0] + "\t" + singles.get(t)[1] + "\t" + singles.get(t)[2]);

            wr.newLine();
        }
        wr.close();
        rd.close();
    }

   

    private void checkLinks(Map<String, String[]> paths, Map<Integer, List<Link>> links, String[] pth, String path_id, int offset, int length) {
        boolean eq;
        List<Link> lk;
        boolean analyzed = false;
        for (Integer it : links.keySet()) {
            if (analyzed) {
                break;
            }
            eq = true;
            lk = links.get(it);
            if (lk.get(0).length == length) {
                for (int i = 0; i < length; i++) {
                    if (!paths.get(lk.get(0).path_id)[lk.get(0).offset + i].equals(pth[i])) {
                        eq = false;
                        break;
                    }
                }
                if (eq) {
                    lk.add(new Link(path_id, offset, length));
                    links.put(it, lk);
                    analyzed = true;
                    break;
                }
            }
        }
        if (!analyzed) {
            lk = new ArrayList();
            lk.add(new Link(path_id, offset, length));
            links.put(links.size() + 1, lk);
        }
    }

    private void writeLinks(Map<String, String[]> paths, Map<Integer, List<Link>> links, String outf) throws IOException {
        List<Link> lk;
        BufferedWriter wr = new BufferedWriter(new FileWriter(outf));
        String line;
        for (Integer it : links.keySet()) {
            lk = links.get(it);
            line = lk.size() + "\t" + lk.get(0).length + "\t";
            for (Link lnk : lk) {
                line = line + lnk.path_id + ";";
            }
            line = line + "\t" + paths.get(lk.get(0).path_id)[1] + "\t";
            line = line + paths.get(lk.get(0).path_id)[paths.get(lk.get(0).path_id).length - 1];
            for (int i = lk.get(0).offset; i < lk.get(0).length + lk.get(0).offset; i++) {
                line = line + "\t" + paths.get(lk.get(0).path_id)[i];
            }
            wr.write(line);
            wr.newLine();
        }
        wr.close();
    }

    

    private void get_connectivity(String infile, String outfile, Network all) throws FileNotFoundException, IOException {
        BufferedReader rd = new BufferedReader(new FileReader(infile));
        BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
        String s;
        String[] ss;
        int len, cn;

        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            len = Integer.parseInt(ss[0]);
            for (int i = 0; i < len; i++) {

            }
        }
    }

    
    private class I2S2 {

        int i1;
        int i2;
        String s1;
        String s2;

        I2S2(int i1, int i2, String s1, String s2) {
            this.s1 = s1;
            this.s2 = s2;
            this.i1 = i1;
            this.i2 = i2;
        }

        @Override
        public String toString() {

            return this.s1 + "\t" + this.s2 + "\t";
        }

        @Override
        public boolean equals(Object s) {
            final I2S2 ss = (I2S2) s;
            if (this.s1.equals(ss.s1) && this.s1.equals(ss.s1) && this.s1.equals(ss.s1)) {
                return true;
            } else {
                return false;
            }

        }

        @Override
        public int hashCode() {
            int hash = 5;
            hash = 83 * hash + this.i1;
            hash = 83 * hash + this.i2;
            hash = 83 * hash + (this.s1 != null ? this.s1.hashCode() : 0);
            hash = 83 * hash + (this.s2 != null ? this.s2.hashCode() : 0);
            return hash;
        }

    }

    private class I1S1 {

        int i1;
        String s1;

        I1S1(int i1, String s1) {
            this.s1 = s1;

            this.i1 = i1;
        }

        @Override
        public String toString() {

            return this.s1 + "\t" + this.i1 + "\t";
        }

        @Override
        public boolean equals(Object s) {
            final I1S1 ss = (I1S1) s;
            if (this.s1.equals(ss.s1) && this.i1 == ss.i1) {
                return true;
            } else {
                return false;
            }

        }

        @Override
        public int hashCode() {
            int hash = 5;
            hash = 83 * hash + this.i1;

            hash = 83 * hash + (this.s1 != null ? this.s1.hashCode() : 0);
            return hash;
        }

    }

    private class S3 {

        String s1;
        String s2;
        String s3;

        S3(String s1, String s2, String s3) {
            this.s1 = s1;
            this.s2 = s2;
            this.s3 = s3;
        }

        @Override
        public String toString() {

            return this.s1 + "\t" + this.s2 + "\t" + this.s3;
        }

        @Override
        public boolean equals(Object s) {
            final S3 ss = (S3) s;
            if (this.s1.equals(ss.s1) && this.s1.equals(ss.s1) && this.s1.equals(ss.s1)) {
                return true;
            } else {
                return false;
            }

        }

        @Override
        public int hashCode() {
            int hash = 3;
            hash = 17 * hash + (this.s1 != null ? this.s1.hashCode() : 0);
            hash = 17 * hash + (this.s2 != null ? this.s2.hashCode() : 0);
            hash = 17 * hash + (this.s3 != null ? this.s3.hashCode() : 0);
            return hash;
        }

    }
}
