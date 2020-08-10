package masterpath.masterpath;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * PathwayManager class contains methods to deal with pathways
 *
 * @author Natalia Rubanova
 */
public class PathwayManager {

    
    
     /**
     * Find the shortest pathways
     *
     * @param foundf Input file
     * @param outname Output file
     * @param all Network
     * @param hg Hit genes list
     * @param fpl Final players list
     * @param d_ppi Length gap for protein-protein interactions
     * @param d_tf Length gap for transcriptional interactions
     * @param d_mirna Length gap for miRNA-mRNA interactions
     * @param d_kegg Length gap for metabolic interactions
     */
    public void find_the_shortest_paths(String foundf, String outname, Network all, Map<String, String> hg, Map<String, String> fpl, int d_ppi, int d_tf, int d_mirna, int d_kegg) throws FileNotFoundException, IOException {
        System.out.println("+++++++++++++find_the_shortest_paths++++++++++++");
        BufferedReader rd = new BufferedReader(new FileReader(foundf));
        BufferedWriter wr = new BufferedWriter(new FileWriter(foundf + outname));
        Map<String, Map<String, List<Path_strength>>> hps = new HashMap(); //hg - path- strength
        List<Path_strength> list_pstr;
        List<String> fplayers_names;
        Map<String, List<Path_strength>> map_pstr;
        Map<String, Map<String, int[]>> map_min_length = new HashMap();
        int[] arr_min_len;
        // 0 - ppi; 1- tf; 2- minra; 3 -kegg;
        String s, type, tmp, hgene, fplayer;
        String[] ss;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            tmp = "";
            if (s.contains("KEGG")) {
                type = "kegg";
            } else {
                if (s.contains("transmir") && (s.contains("MIRT") || s.contains("sl:")) || s.contains("_mr")) {
                    type = "mirna";
                } else {
                    for (int i = 3; i < ss.length; i++) {
                        tmp = tmp + "\t" + ss[i];
                    }
                    if (tmp.contains("man") || tmp.contains("tfacts:") || tmp.contains("XN000") || tmp.contains("_tf")) {
                        type = "tf";
                    } else {
                        type = "ppi";
                    }
                }
            }
            hgene = ss[1];
            fplayer = ss[ss.length - 1];
            if (hps.containsKey(hgene)) {
                if (hps.get(hgene).containsKey(fplayer)) {
                    /* list_pstr = hps.get(hgene).get(fplayer);
                     list_pstr.add(new Path_strength(s, new Strength(type, ss.length - 3)));
                     hps.get(hgene).put(fplayer, list_pstr); */

                    if (type.equals("ppi") & ss.length - 3 <= map_min_length.get(hgene).get(fplayer)[0] + d_ppi) {
                        list_pstr = hps.get(hgene).get(fplayer);
                        list_pstr.add(new Path_strength(s, new Strength(type, ss.length - 3)));
                        hps.get(hgene).put(fplayer, list_pstr);
                    }
                    if (type.equals("tf") & ss.length - 3 <= map_min_length.get(hgene).get(fplayer)[1] + d_tf) {
                        list_pstr = hps.get(hgene).get(fplayer);
                        list_pstr.add(new Path_strength(s, new Strength(type, ss.length - 3)));
                        hps.get(hgene).put(fplayer, list_pstr);
                    }
                    if (type.equals("mirna") & ss.length - 3 <= map_min_length.get(hgene).get(fplayer)[2] + d_mirna) {
                        list_pstr = hps.get(hgene).get(fplayer);
                        list_pstr.add(new Path_strength(s, new Strength(type, ss.length - 3)));
                        hps.get(hgene).put(fplayer, list_pstr);
                    }
                    if (type.equals("kegg") & ss.length - 3 <= map_min_length.get(hgene).get(fplayer)[3] + d_kegg) {
                        list_pstr = hps.get(hgene).get(fplayer);
                        list_pstr.add(new Path_strength(s, new Strength(type, ss.length - 3)));
                        hps.get(hgene).put(fplayer, list_pstr);
                    }

                    if (type.equals("ppi") && ss.length - 3 < map_min_length.get(hgene).get(fplayer)[0]) {
                        map_min_length.get(hgene).get(fplayer)[0] = ss.length - 3;
                        /* tmp_min_len = min_length.get(hgene);
                         tmp_min_len[0] = ss.length - 3;
                         min_length.put(hgene, tmp_min_len); */
                    }
                    if (type.equals("tf") && ss.length - 3 < map_min_length.get(hgene).get(fplayer)[1]) {
                        map_min_length.get(hgene).get(fplayer)[1] = ss.length - 3;
                        /* tmp_min_len = min_length.get(hgene);
                         tmp_min_len[1] = ss.length - 3;
                         min_length.put(hgene, tmp_min_len);*/
                    }

                    if (type.equals("mirna") && ss.length - 3 < map_min_length.get(hgene).get(fplayer)[2]) {
                        map_min_length.get(hgene).get(fplayer)[2] = ss.length - 3;
                        /*  tmp_min_len = min_length.get(hgene);
                         tmp_min_len[2] = ss.length - 3;
                         min_length.put(hgene, tmp_min_len);*/
                    }
                    if (type.equals("kegg") && ss.length - 3 < map_min_length.get(hgene).get(fplayer)[3]) {
                        map_min_length.get(hgene).get(fplayer)[3] = ss.length - 3;
                        /*  tmp_min_len = min_length.get(hgene);
                         tmp_min_len[3] = ss.length - 3;
                         min_length.put(hgene, tmp_min_len);*/
                    }
                } else {
                    list_pstr = new ArrayList();
                    list_pstr.add(new Path_strength(s, new Strength(type, ss.length - 3)));
                    hps.get(hgene).put(fplayer, list_pstr);

                    arr_min_len = new int[]{Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE};

                    if (type.equals("ppi")) {
                        arr_min_len[0] = ss.length - 3;
                        map_min_length.get(hgene).put(fplayer, arr_min_len);
                    }
                    if (type.equals("tf")) {
                        arr_min_len[1] = ss.length - 3;
                        map_min_length.get(hgene).put(fplayer, arr_min_len);
                    }
                    if (type.equals("mirna")) {
                        arr_min_len[2] = ss.length - 3;
                        map_min_length.get(hgene).put(fplayer, arr_min_len);
                    }
                    if (type.equals("kegg")) {
                        arr_min_len[3] = ss.length - 3;
                        map_min_length.get(hgene).put(fplayer, arr_min_len);
                    }
                }
            } else {
                list_pstr = new ArrayList();
                list_pstr.add(new Path_strength(s, new Strength(type, ss.length - 3)));
                hps.put(hgene, new HashMap());
                hps.get(hgene).put(fplayer, list_pstr);

                arr_min_len = new int[]{Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE};
                if (type.equals("ppi")) {
                    arr_min_len[0] = ss.length - 3;
                    map_min_length.put(hgene, new HashMap());
                    map_min_length.get(hgene).put(fplayer, arr_min_len);
                }
                if (type.equals("tf")) {
                    arr_min_len[1] = ss.length - 3;
                    map_min_length.put(hgene, new HashMap());
                    map_min_length.get(hgene).put(fplayer, arr_min_len);
                }
                if (type.equals("mirna")) {
                    arr_min_len[2] = ss.length - 3;
                    map_min_length.put(hgene, new HashMap());
                    map_min_length.get(hgene).put(fplayer, arr_min_len);
                }
                if (type.equals("kegg")) {
                    arr_min_len[3] = ss.length - 3;
                    map_min_length.put(hgene, new HashMap());
                    map_min_length.get(hgene).put(fplayer, arr_min_len);
                }
            }
        }

        for (String hgn : hps.keySet()) {
            fplayers_names = new ArrayList<String>(hps.get(hgn).keySet());
            Collections.sort(fplayers_names);
            System.out.println(hgn);
            for (String fplr : fplayers_names) {
                System.out.print("\t" + fplr + "\t");
                for (int j = 0; j < map_min_length.get(hgn).get(fplr).length; j++) {
                    if (map_min_length.get(hgn).get(fplr)[j] < Integer.MAX_VALUE) {
                        System.out.print(map_min_length.get(hgn).get(fplr)[j] + "\t"); 
                    } else {
                        System.out.print("NA" + "\t");
                    }
                }
                System.out.print("\n");

                // System.out.println(t + "\t" + min_length.get(t)[0] + "\t" + min_length.get(t)[1] + "\t" + min_length.get(t)[2] + "\t" + min_length.get(t)[3]);
                for (Path_strength ps : hps.get(hgn).get(fplr)) {
                    if (ps.strength.type.equals("ppi") & ps.strength.len <= map_min_length.get(hgn).get(fplr)[0] + d_ppi) {
                        wr.write(ps.path + "\n");
                    }
                    if (ps.strength.type.equals("tf") & ps.strength.len <= map_min_length.get(hgn).get(fplr)[1] + d_tf) {
                        wr.write(ps.path + "\n");
                    }
                    if (ps.strength.type.equals("mirna") & ps.strength.len <= map_min_length.get(hgn).get(fplr)[2] + d_mirna) {
                        wr.write(ps.path + "\n");
                    }
                    if (ps.strength.type.equals("kegg") & ps.strength.len <= map_min_length.get(hgn).get(fplr)[3] + d_kegg) {
                        wr.write(ps.path + "\n");
                    }
                }
            }
        }
        rd.close();
        wr.close();
    }
    
    /**
     * Find the strongest pathways
     *
     * @param foundf Input file
     * @param outname Output file     
     */

    public void find_the_strongest_pathways(String foundf, String outfile) throws IOException {
        BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
        BufferedReader rd = new BufferedReader(new FileReader(foundf));

        String s;
        String[] ss;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            if (Float.parseFloat(ss[1]) <= 0.002 && Float.parseFloat(ss[0]) >= 2) {
                wr.write(ss[3]);
                for (int i = 4; i < ss.length; i++) {
                    wr.write("\t" + ss[i] );
                }
                wr.write( "\n");
            }
        }
        rd.close();
        wr.close();

    }
    
    
    
    /**
     * Find and rank miRNAs inside pathways
     *
     * @param foundf File with common paths
     * @param resf Output file
     * @param mirtarbase Link to mirTarBase network
     * @param hg Hit genes list
     * @param fpl Final implementer list
     * @param length Maximum length
     * @param min
     * @param mask
     */
    public void find_miRNAs_on_pathways(String foundf, String resf, Network mirtarbase, Map<String, String> hg, Map<String, String> fpl, int length, int min, String mask) throws FileNotFoundException, IOException {
        System.out.println("+++++++++++++calcmiRNAs++++++++++++");
        BufferedReader rd = new BufferedReader(new FileReader(foundf));
        //BufferedWriter wr= new BufferedWriter (new FileWriter(resf));
        String s;//, mi;
        String[] ss;
        int count;
        Map<String, Integer> id_count = new HashMap();
        Map<String, int[]> id_count_by_3 = new HashMap();        //mirtarbase id - data
        Map<String, List<Map<String, Integer>>> id_list = new HashMap();

        Map<String, Integer> mi_count = new HashMap(); //mirna - number of occurrence
        Map<String, int[]> mi_count_by_3 = new HashMap();        //mirna - data
        Map<String, List<Map<String, Integer>>> mi_list = new HashMap();

        int[] ln;
        List<Map<String, Integer>> ln_m;
        Map<String, Integer> ln_m_tmp;
        String beg, end;
        List<String> all_mi;
        int d;
        int start_pos;
        String mi;
        while ((s = rd.readLine()) != null) {
            if (s.contains("MIRT")) {
                all_mi = new ArrayList();
                ss = s.split("\t");
                if (mask.equals("countonlymiddle")) {
                    start_pos = 3;
                } else {
                    start_pos = 2;
                }
                for (int i = start_pos; i < ss.length; i++) {
                    if (ss[i].contains("MIRT")) {
                        all_mi.add(ss[i]);
                    }
                }
                for (String mir_int : all_mi) {
                    mi = mirtarbase.interactions.get(mir_int).int1.id;
                    if (id_count.containsKey(mi)) {
                        ////////////////////////////////
                        count = id_count.get(mi);
                        id_count.put(mi, count + 1);
                        ////////////////////////////////
                        ss = s.split("\t");
                        ln = id_count_by_3.get(mi);
                        for (int i = 0; i < length - min; i++) {
                            if (ss.length == i + 5) {
                                ln[i]++;
                            }
                        }
                        id_count_by_3.put(mi, ln);

                        ////////////////////////////////
                        //*1
                        beg = hg.get(ss[1]);
                        end = fpl.get(ss[ss.length - 1]);
                        ln_m = id_list.get(mi);
                        for (int i = 0; i < length - min; i++) {
                            if (ss.length == i + 6) {
                                ln_m_tmp = ln_m.get(i);
                                if (ln_m.get(i).containsKey(beg + "-" + end)) {
                                    d = ln_m_tmp.get(beg + "-" + end) + 1;
                                    ln_m_tmp.put(beg + "-" + end, d);
                                } else {
                                    ln_m_tmp.put(beg + "-" + end, 1);
                                }
                                ln_m.set(i, ln_m_tmp);
                                break;
                            }
                        }
                        id_list.put(mi, ln_m);
                        ////////////////////////////////
                    } else {
                        ////////////////////////////////
                        id_count.put(mi, 1);
                        //////////////////////////////
                        ln = new int[length - min];
                        /*  for (int i = 0; i < length - min; i++) {
                         ln[i] = 0;
                         }
                         ss = s.split("\t"); */
                        for (int i = 0; i < length - min; i++) {
                            ln[i] = 0;
                            if (ss.length == i + 5) {
                                ln[i] = 1;
                            }
                        }
                        id_count_by_3.put(mi, ln);
                        //////////////////////////////////*2
                        beg = hg.get(ss[1]);  // bez + beg\end
                        end = fpl.get(ss[ss.length - 1]);
                        ln_m = new ArrayList();
                        ln_m_tmp = new HashMap();
                        for (int i = 0; i < length - min; i++) {
                            ln_m.add(new HashMap());
                            if (ss.length == i + 6) {
                                ln_m_tmp.put(beg + "-" + end, 1);
                                ln_m.set(i, ln_m_tmp);
                            }
                        }
                        id_list.put(mi, ln_m);
                    }
                }
            }
        }

        count = 0;
        for (String t : id_count.keySet()) {
            count = count + id_count.get(t);
            System.out.print(t + "\t" + id_count.get(t) + "\t");
            for (int i = 0; i < length - min; i++) {
                System.out.print(id_count_by_3.get(t)[i] + "\t");
            }
            for (int i = 0; i < length - min; i++) {
                System.out.print(i + length - min - 1 + ":\t");
                for (String tt : id_list.get(t).get(i).keySet()) {
                    System.out.print(tt + " x " + id_list.get(t).get(i).get(tt) + "\t");
                }
            }
            System.out.println();
        }
        System.out.println("count " + count);
        rd.close();
    }

    /**
     * Find and rank miRNAs inside paths
     *
     * @param foundf Common paths file
     * @param resf Output file
     * @param mirtarbase Link to mirTarBAse network
     * @param transmir Link to transMir network
     * @param hg Hit genes list
     * @param fpl Final players List
     * @param length Maximum length
     * @param min_occ Minimum occurrence
     * @param mask Mask
     */
    public void find_miRNAs_on_paths(String foundf, String resf, Network mirtarbase, Network transmir, Map<String, String> hg, Map<String, String> fpl, int length, int min_occ, String mask) throws FileNotFoundException, IOException {
        System.out.println("+++++++++++++calcmiRNAs_overrepr++++++++++++");
        BufferedReader rd = new BufferedReader(new FileReader(foundf));
        //BufferedWriter wr= new BufferedWriter (new FileWriter(resf));
        String s;
        String[] ss;
        Map<String, List<String[]>> mirnas = new HashMap();
        TreeMap<String, Integer> mirnas_maxocc = new TreeMap();
        List<String_Int> mirnas_maxocc_list = new ArrayList();
        List<String[]> list;
        String mi, t;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            if (Integer.parseInt(ss[0]) >= min_occ) {
                for (int i = 5; i < ss.length - Integer.parseInt(ss[1]); i++) {
                    if (ss[i].contains("MIRT")) {
                        mi = mirtarbase.interactions.get(ss[i]).int1.id;
                        if (mirnas.containsKey(mi)) {
                            list = mirnas.get(mi);
                            list.add(ss);
                            mirnas.put(mi, list);
                            if (Integer.parseInt(ss[0]) > mirnas_maxocc.get(mi)) {
                                mirnas_maxocc.put(mi, Integer.parseInt(ss[0]));
                            }
                        } else {
                            list = new ArrayList();
                            list.add(ss);
                            mirnas.put(mi, list);
                            mirnas_maxocc.put(mi, Integer.parseInt(ss[0]));
                        }
                    }
                    if (ss[i].contains("transmir")) {
                        mi = transmir.interactions.get(ss[i]).int2.id;
                        if (mirnas.containsKey(mi)) {
                            list = mirnas.get(mi);
                            list.add(ss);
                            mirnas.put(mi, list);
                            if (Integer.parseInt(ss[0]) > mirnas_maxocc.get(mi)) {
                                mirnas_maxocc.put(mi, Integer.parseInt(ss[0]));
                            }
                        } else {
                            list = new ArrayList();
                            list.add(ss);
                            mirnas.put(mi, list);
                            mirnas_maxocc.put(mi, Integer.parseInt(ss[0]));
                        }
                    }
                }
            }
        }

        for (String tt : mirnas_maxocc.keySet()) {
            mirnas_maxocc_list.add(new String_Int(tt, mirnas_maxocc.get(tt)));
        }
        Collections.sort(mirnas_maxocc_list, new MaxOccComparator());

        for (String_Int strint : mirnas_maxocc_list) {
            t = strint.str;
            System.out.println(t + "\t" + strint.it);
            Collections.sort(mirnas.get(t), new ListComparator());
            for (String[] get : mirnas.get(t)) {
                System.out.print("--");
                for (String get1 : get) {
                    System.out.print("\t" + get1);
                }
                System.out.println();
            }
        }
        rd.close();
    }

    /**
     * Create a file for Cytoscape software from the shortest paths
     *
     * @param foundf File with paths
     * @param all Network
     * @param hg Hit genes list
     * @param fpl Final players list
     * @throws FileNotFoundException
     * @throws IOException
     */
    public void create_cyto_for_paths(String foundf, Network all, Map<String, String> hg, Map<String, String> fpl) throws FileNotFoundException, IOException {
        BufferedReader rd = new BufferedReader(new FileReader(foundf));
        BufferedWriter wr = new BufferedWriter(new FileWriter(foundf + "_cyto"));
        String s, fp;
        String[] ss;
        int count = 0;
        int id_path = 0;
        int ss_len;
        //String[] fpls = {"pHGNC:7611", "pHGNC:4223", "pHGNC:7566", "pHGNC:5466"};
        Map<String, Integer> fpl_count = new HashMap();
        /* fpl_count.put("pHGNC:7611", 0);
         fpl_count.put("pHGNC:4223", 0);
         fpl_count.put("pHGNC:7566", 0);
         fpl_count.put("pHGNC:5466", 0); */

        for (String t : fpl.keySet()) {
            fpl_count.put(t, 0);
        }
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            ss_len = ss.length;
            fp = ss[4];
            count = fpl_count.get(fp);
            fpl_count.put(fp, count + 1);
            for (int i = 5; i < ss.length; i++) {
                wr.write(fp + "_" + fpl_count.get(fp) + "_" + i + "\t");
                wr.write(ss[2] + "\t");
                wr.write(ss[3] + "\t");
                wr.write(ss[0] + "\t");
                wr.write(ss[1] + "\t");
                wr.write(ss[i] + "\t");
                try {
                    wr.write(all.interactions.get(ss[i]).int1.id + "\t");
                } catch (NullPointerException e) {
                    System.out.println(ss[i]);
                    System.out.println(all.interactions.get(ss[i]).sourcedbentry.get(0));
                    break;
                }
                wr.write(all.interactions.get(ss[i]).int2.id + "\t");
                wr.write(fp + "_" + fpl_count.get(fp) + "\t");
                wr.newLine();
            }
        }
        rd.close();
        wr.close();
    }

    /**
     * Creates a file for Cytoscape software for a list of pathways by pathways
     * ids
     *
     * @param path_file File with pathways
     * @param filename File with list of pathways by id
     * @param out_dir Output folder
     * @param all Network
     * @param hg Hit genes list
     * @param fpl Final players list
     */
    public void create_cyto_for_pathways(String path_file, String filename, String out_dir, Network all, Map<String, String> hg, Map<String, String> fpl) throws FileNotFoundException, IOException {

        BufferedReader rd = new BufferedReader(new FileReader(path_file));
        BufferedWriter wr;
        Map<String, Integer> fpl_count = new HashMap();
        Map<String, String> paths = new HashMap();
        String s, fp, pth;
        String[] ss, ids, pth_ss;
        int ss_len, pth_len, count = 0, id_path = 0;

        for (String t : fpl.keySet()) {
            fpl_count.put(t, 0);
            // System.out.println(t);
        }
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            paths.put(ss[0], s);
        }
        rd.close();
        rd = new BufferedReader(new FileReader(filename));
        while ((s = rd.readLine()) != null) {

            ss = s.split("\t");
            wr = new BufferedWriter(new FileWriter(out_dir + ss[0] + "_cyto"));
            ids = ss[1].split(";");
            for (int j = 0; j < ids.length; j++) {
                pth = paths.get(ids[j]);
                pth_ss = pth.split("\t");
                pth_len = pth_ss.length;
                fp = pth_ss[pth_len - 1];
                // System.out.println("------- " + fp);
                count = fpl_count.get(fp);
                fpl_count.put(fp, count + 1);
                for (int i = 2; i < pth_len - 1; i++) {
                    wr.write(fp + "_" + fpl_count.get(fp) + "_" + i + "\t");
                    wr.write(pth_ss[2] + "\t");
                    wr.write(pth_ss[3] + "\t");
                    wr.write(pth_ss[0] + "\t");
                    wr.write(pth_ss[1] + "\t");
                    wr.write(pth_ss[i] + "\t");
                    try {
                        wr.write(all.interactions.get(pth_ss[i]).int1.id + "\t");
                    } catch (NullPointerException e) {
                        System.out.println(pth_ss[i]);
                        System.out.println(ids[j]);
                        System.out.println(all.interactions.get(pth_ss[i]).sourcedbentry.get(0));
                        break;
                    }
                    wr.write(all.interactions.get(pth_ss[i]).int2.id + "\t");
                    wr.write(fp + "_" + fpl_count.get(fp) + "\t");
                    wr.newLine();
                }
            }
            wr.close();
        }
        rd.close();
    }
    
        /**
     * Creates a file for Cytoscape software for a list of pathways 
     * 
     *
     * @param path_file File with pathways
     * @param filename File with list of pathways by id
     * @param out_dir Output folder
     * @param all Network
     * @param hg Hit genes list
     * @param fpl Final players list
     */
    
    public void create_cyto_for_list_of_pathways(String path_file, String filename, String out_dir, Network all, Map<String, String[]> hugo_by_id) throws FileNotFoundException, IOException {

        BufferedReader rd ;
        BufferedWriter wr_nodes = new BufferedWriter(new FileWriter(filename + "_nodes_cyto"));
        BufferedWriter wr_interactions = new BufferedWriter(new FileWriter(filename + "_intercations_cyto"));

        String s;
        String[] ss;



  
        rd = new BufferedReader(new FileReader(filename));
        ArrayList<String> pathsids = new ArrayList();
        
        while ((s = rd.readLine()) != null) {
            if (!pathsids.contains(s)){
                pathsids.add(s);
            }
        }
        
        rd = new BufferedReader(new FileReader(path_file));
        
        String id, intr, intr1, intr2, type, dir, entry, symb1, symb2,type1, type2, entry1, entry2;
        ArrayList<String> all_intr = new ArrayList();
        ArrayList<String> all_nodes = new ArrayList();
       
        while ((s = rd.readLine()) != null) {     
            ss = s.split("\t");
            id = ss[0];            
            if (pathsids.contains(id)) {
                for ( int j = 2; j < ss.length -1 ; j++){
                    intr = ss[j];
                    try {
                        intr1 = all.interactions.get(intr).int1.id;                      
                    } catch (NullPointerException e) {
                        intr1 = all.interactions.get(intr).sourcedbentry.get(0);
                    }
                    try {
                        intr2 = all.interactions.get(intr).int2.id;                      
                    } catch (NullPointerException e) {
                        intr2 = all.interactions.get(intr).sourcedbentry.get(0);
                    }
                    dir = all.interactions.get(intr).dir;
                    type = all.interactions.get(intr).type;
                    
                    entry = intr1 + "\t" + intr2 + "\t" + type + "\t" + dir;
                    if (!all_intr.contains(entry)){
                        all_intr.add(entry);
                    }
                    
                    if (intr1.startsWith("pHGNC")){
                        symb1 = "p" + hugo_by_id.get(intr1.substring(1))[1];
                    } else if (intr1.startsWith("HGNC")) {
                        symb1 = hugo_by_id.get(intr1)[1];
                    } else {
                        symb1 = intr1;
                    }
                    if (intr2.startsWith("pHGNC")){
                        symb2 = "p" + hugo_by_id.get(intr2.substring(1))[1];
                    } else if (intr1.startsWith("HGNC")) {
                        symb2 = hugo_by_id.get(intr2)[1];
                    } else {
                        symb2 = intr2;
                    }
                    
                    type1 = all.nodes.get(intr1).type;
                    
                    type2 = all.nodes.get(intr2).type;
                    entry1 = intr1 + "\t" + symb1 + "\t" + type1;
                    entry2 = intr2 + "\t" + symb2 + "\t" + type2;
                    if (!all_nodes.contains(entry1)){
                        all_nodes.add(entry1);
                    }
                    if (!all_nodes.contains(entry2)){
                        all_nodes.add(entry2);
                    }                    
                }
            } 
        }
        
        for (String ent : all_intr){
            wr_interactions.write(ent + "\n");
        }
        for (String ent : all_nodes){
            wr_nodes.write(ent + "\n");
        }
        wr_interactions.close();
        wr_nodes.close();
        rd.close();
    }
        
        
    private void create_cyto_for_comparisson(String com, String p1, String p1_single_max, String p2, String p2_single_max, String out, Network all) throws FileNotFoundException, IOException {
        //, Map<String, String> hg, Map<String, String> fpl
        BufferedReader rd;
        BufferedWriter wr = new BufferedWriter(new FileWriter(out));
        String s;
        String[] ss, tt;
        int off, len, j;
        Set<String> p_id1 = new HashSet();
        Set<String> p_id2 = new HashSet();
        Set<String> pth = new HashSet();
        //Set <String> pth2= new HashSet();
        Map<String, String> comp = new HashMap();
        Map<String, String[]> orp_pth1 = new HashMap();
        Map<String, String[]> orp_pth2 = new HashMap();
        Map<String, String> cyto = new HashMap();
        rd = new BufferedReader(new FileReader(com));
        while ((s = rd.readLine()) != null) {
            off = 0;
            ss = s.split("\t");
            for (String s1 : ss) {
                if ("".equals(s1)) {
                    off++;
                }
            }
            len = (ss.length - 1 - off) / 2;
            j = 0;
            for (int i = 1 + 5; i < len + 1; i++) {
                j++;
                //cyto.put(ss[0] + "_" + j, all.interactions.get(ss[i]).int1.id + "\t" + all.interactions.get(ss[i]).int2.id + "\t" + ss[1] + "\t" + ss[1 + len + off] + "\t" + ss[2]);
                if (comp.containsKey(ss[i]) && Integer.parseInt(comp.get(ss[i]).split("\t")[1]) + Integer.parseInt(comp.get(ss[i]).split("\t")[2]) < Integer.parseInt(ss[1] + Integer.parseInt(ss[1 + len + off]))) {
                    comp.put(ss[i], ss[0] + "_" + j + "\t" + ss[1] + "\t" + ss[1 + len + off] + "\t" + ss[2] + "\t" + "sh");
                } else {
                    comp.put(ss[i], ss[0] + "_" + j + "\t" + ss[1] + "\t" + ss[1 + len + off] + "\t" + ss[2] + "\t" + "sh");
                }
            }
            tt = ss[3].split(";");
            for (String t : tt) {
                if (!"".equals(t)) {
                    p_id1.add(t);
                }
            }
            tt = ss[1 + len + off + 2].split(";");
            for (String t : tt) {
                if (!"".equals(t)) {
                    p_id2.add(t);
                }
            }
        }
        rd = new BufferedReader(new FileReader(p1_single_max));
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            orp_pth1.put(ss[0], ss);
        }
        rd = new BufferedReader(new FileReader(p2_single_max));
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            orp_pth2.put(ss[0], ss);
        }
        rd = new BufferedReader(new FileReader(p1));
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            if (p_id1.contains(ss[0])) {
                pth.add(s);
            }
        }
        rd = new BufferedReader(new FileReader(p2));
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            if (p_id2.contains(ss[0])) {
                pth.add(s);
            }
        }

        for (String t : pth) {
            ss = t.split("\t");
            j = 0;
            for (int i = 2; i < ss.length - 1; i++) {
                j++;
                if (!comp.containsKey(ss[i])) {
                    //System.out.println(ss[i]);
                    if (orp_pth1.containsKey(ss[i])) {
                        if (orp_pth2.containsKey(ss[i])) {
                            comp.put(ss[i], ss[0] + "_" + j + "\t" + orp_pth1.get(ss[i])[1] + "\t" + orp_pth2.get(ss[i])[1] + "\t" + orp_pth1.get(ss[i])[2] + "\t" + "12");
                        } else {
                            comp.put(ss[i], ss[0] + "_" + j + "\t" + orp_pth1.get(ss[i])[1] + "\t" + 0 + "\t" + orp_pth1.get(ss[i])[2] + "\t" + "1");
                        }
                    } else {
                        comp.put(ss[i], ss[1] + "_" + j + "\t" + 0 + "\t" + orp_pth2.get(ss[i])[1] + "\t" + orp_pth2.get(ss[i])[2] + "\t" + "2");
                    }
                }
            }
        }

        for (String t : comp.keySet()) {
            ss = comp.get(t).split("\t");
            wr.write(ss[0] + "\t" + all.interactions.get(t).int1.id + "\t" + all.interactions.get(t).int2.id + "\t" + ss[1] + "\t" + ss[2] + "\t" + ss[3] + "\t" + ss[4]);
            wr.newLine();
        }

        rd.close();
        wr.close();
    }

    /**
     * Filter pathways by length
     *
     * @param inf Input file
     * @param outf Output file
     * @param min Minimum length
     * @param max Maximum length
     */
    public void filterPathways_by_length(String inf, String outf, int min, int max) throws IOException {
        BufferedReader rd = new BufferedReader(new FileReader(inf));
        BufferedWriter wr = new BufferedWriter(new FileWriter(inf + outf));
        String s;
        String[] ss;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            if (ss.length >= min + 3 && ss.length <= max + 3) {
                wr.write(s);
                wr.newLine();
            }
        }
        rd.close();
        wr.close();
    }
        /**
     * Filter pathways by length
     *
     * @param inf Input file
     * @param outf Output file
     * @param min Minimum length
     * @param max Maximum length
     */
    public void filterPathways_by_centrality_score(String inf, String outf, int min, int max) throws IOException {
        BufferedReader rd = new BufferedReader(new FileReader(inf));
        BufferedWriter wr = new BufferedWriter(new FileWriter(outf));
        String s;
        String[] ss;
        int score;
        int length;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            score =Integer.parseInt(ss[6]);
            length =Integer.parseInt(ss[0]);
            if ((length>2 && score >= min )) {
                //(length ==1 && score>=10 ) || 
                wr.write(s);
                wr.newLine();
            }
        }
        rd.close();
        wr.close();
    }
    
    /**
     * Filter pathways by node names \ type
     *
     * @param foundf Input file
     * @param outname Output file
     * @param all Network
     * @param hg Hit genes list
     * @param fpl Final players list
     * @param mask Mask
     * @param gene Node name
     * @param type Node type
     */
    public void filterPathways(String foundf, String outname, Network all, Map<String, String> hg, Map<String, String> fpl, String mask, String gene, String type) throws FileNotFoundException, IOException {
        System.out.println("+++++++++++++filter pathways++++++++++++");
        BufferedReader rd = new BufferedReader(new FileReader(foundf));
        BufferedWriter wr = new BufferedWriter(new FileWriter(outname));
        String s, gene_name;
        String[] ss;
        List<String> checked_mirnas;

        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            if (mask.equals("mir_middle") && s.contains("transmir")) { //mirna in the middle
                wr.write(s);
                wr.newLine();
            }
            if (mask.equals("mir_all") && s.contains("MIRT")) { //mirna on the pathway including in the begining
                wr.write(s);
                wr.newLine();
            }
            if (mask.equals("hg") && ss[1].equals(gene) && s.contains("MIRT")) { // starts with specific hg, mirna on the pathway
                wr.write(s);
                wr.newLine();
            }
            if (mask.equals("mirna") && s.contains("MIRT")) { //mirna on the pathway , pathway contains specific mirna
                ss = s.split("\t");
                checked_mirnas = new ArrayList();
                for (int i = 2; i < ss.length - 1; i++) {
                    gene_name = all.interactions.get(ss[i]).int1.id;
                    if (ss[i].contains("MIRT") && gene_name.equals(gene)) {
                        if (!checked_mirnas.contains(gene_name)) {
                            wr.write(s + "\n");
                            checked_mirnas.add(gene_name);
                        }
                    }
                }
            }

        }
        rd.close();
        wr.close();
    }

    /**
     * Filter paths on condition ppi [3,4] & miRNA [6/5] or ppi [3,4] or miRNA
     * [6/5]
     *
     * @param foundf File with pathways
     * @param outname Output file
     * @param all Network
     * @param hg Hit genes
     * @param fpl Final implementers
     * @param gene Gene name
     * @param l1 min length for ppi path
     * @param l2 max length for ppi path
     * @param l3 min length for miRNA path
     * @param l4 max length for miRNA path
     */
    public void filterPathways_2(String foundf, String outname, Network all, Map<String, String> hg, Map<String, String> fpl, String gene, int l1, int l2, int l3, int l4) throws FileNotFoundException, IOException {
        //filter on condition ppi [3,4] & miRNA [6/5] or ppi [3,4] or miRNA [6/5]
        System.out.println("+++++++++++++filter pathways_2++++++++++++");
        BufferedReader rd = new BufferedReader(new FileReader(foundf));
        BufferedWriter wr = new BufferedWriter(new FileWriter(foundf + outname));
        String s;
        String[] ss;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            if (s.contains("transmir") && ss.length <= l4 + 3 && ss.length >= l3 + 3) {
                for (int i = 2; i < ss.length - 1; i++) {
                    if (ss[i].contains("MIRT") && all.interactions.get(ss[i]).int1.id.equals(gene)) {
                        wr.write(s + "\n");
                        break;
                    }
                }
            } else {
                if (ss[1].equals(gene) && ss.length <= l2 + 3 && ss.length >= l1 + 3) {
                    wr.write(s);
                    wr.newLine();
                }
            }
        }
        rd.close();
        wr.close();
    }

   

    /**
     * Calculate pathways connectivity
     *
     * @param all Network
     * @param infile Input file
     * @param outfile Output file
     */
    public void add_connectivity_to_pathways(Network all, String infile, String outfile) throws IOException {
        BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
        BufferedReader rd = new BufferedReader(new FileReader(infile));

        int con, con1, con2;
        String s, scon, scon1, scon2;
        String[] ss;

        while ((s = rd.readLine()) != null) {
            con = 0;
            scon = "";
            ss = s.split("\t");
            for (int j = 2; j < ss.length - 1; j++) {
                try {
                    con1 = (all.interactions.get(ss[j]).int1.downnbrs.size() + all.interactions.get(ss[j]).int1.revnbrs.size()) * (all.interactions.get(ss[j]).int1.upnbrs.size() + all.interactions.get(ss[j]).int1.revnbrs.size());
                    con2 = (all.interactions.get(ss[j]).int2.downnbrs.size() + all.interactions.get(ss[j]).int2.revnbrs.size()) * (all.interactions.get(ss[j]).int2.upnbrs.size() + all.interactions.get(ss[j]).int2.revnbrs.size());
                    scon1 = "" + (all.interactions.get(ss[j]).int1.downnbrs.size() + all.interactions.get(ss[j]).int1.revnbrs.size()) + "*" + (all.interactions.get(ss[j]).int1.upnbrs.size() + all.interactions.get(ss[j]).int1.revnbrs.size()) + " " + all.interactions.get(ss[j]).int1.id;
                    scon2 = "" + (all.interactions.get(ss[j]).int2.downnbrs.size() + all.interactions.get(ss[j]).int2.revnbrs.size()) + "*" + (all.interactions.get(ss[j]).int2.upnbrs.size() + all.interactions.get(ss[j]).int2.revnbrs.size()) + " " + all.interactions.get(ss[j]).int2.id;
                } catch (NullPointerException e) {
                    con1 = 1;
                    con2 = 1;
                    scon1 = "*";
                    scon2 = "*";
                }
                if (con1 > con) {
                    con = con1;
                    scon = scon1;
                }
                if (con2 > con) {
                    con = con2;
                    scon = scon2;
                }
            }
            wr.write(s + "\t" + con + "\t" + scon + "\n");
        }

        wr.close();

    }

    /**
     * Filter pathways by connectivity
     *
     * @param in Input file
     * @param out Output file
     * @param out_for_overrreps Output file only for pathway without additional
     * information
     * @param in_min Minimum inward connectivity
     * @param in_max Maximum inward connectivity
     * @param out_min Minimum outward connectivity
     * @param out_max Maximum outward connectivity
     * @throws FileNotFoundException
     * @throws IOException
     */
    public void filter_pathways_by_connectivity(String in, String out, String out_for_overrreps, int in_min, int in_max, int out_min, int out_max) throws FileNotFoundException, IOException {
        BufferedReader rd = new BufferedReader(new FileReader(in));
        BufferedWriter wr = new BufferedWriter(new FileWriter(out));
        BufferedWriter wr_2 = new BufferedWriter(new FileWriter(out_for_overrreps));
        String s;
        String[] ss, tt, rr;
        int in_c, out_c;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            tt = ss[ss.length - 1].split(" ");
            rr = tt[0].split("\\*");
            in_c = Integer.valueOf(rr[0]);
            out_c = Integer.valueOf(rr[1]);
            if ((in_min < in_c) && (in_max > in_c) && (out_min < out_c) && (out_max > out_c)) {
                wr.write(s);
                wr.newLine();
                for (int i = 0; i < ss.length - 3; i++) {
                    wr_2.write(ss[i] + "\t");
                }
                wr_2.write(ss[ss.length - 3] + "\n");
            }
        }
        rd.close();
        wr.close();
        wr_2.close();
    }

    /**
     * Compare two files with paths
     *
     * @param all Network
     * @param hugo_by_id Map with HGNC ids
     * @param f1 File1
     * @param f2 File 2
     * @param out OUtput file
     */
    public void compare_two_paths_files(Network all, Map<String, String[]> hugo_by_id, String f1, String f2, String out) throws FileNotFoundException, IOException {
        BufferedReader rd1 = new BufferedReader(new FileReader(f1));
        BufferedReader rd2 = new BufferedReader(new FileReader(f2));
        BufferedWriter wr = new BufferedWriter(new FileWriter(out));
        String s, path, tmp, n1, n2, names;
        String[] ss, tt;
        Map<String, String> f1_paths = new HashMap();
        Map<String, String> f2_paths = new HashMap();

        int max_len1 = 0, max_len2 = 0;
        int count = 0;
        while ((s = rd1.readLine()) != null) {
            ss = s.split("\t");
            if (Float.parseFloat(ss[ss.length - 1]) <= 0.05 && Integer.parseInt(ss[0]) > 1) {
                path = ss[1];
                for (int i = 2; i <= Integer.parseInt(ss[0]); i++) {
                    path = path + "\t" + ss[i];
                }

                f1_paths.put(path, s);

                if (max_len1 < Integer.parseInt(ss[0])) {
                    max_len1 = Integer.parseInt(ss[0]);
                }
            }
        }
        while ((s = rd2.readLine()) != null) {
            ss = s.split("\t");
            path = ss[1];
            if (Float.parseFloat(ss[ss.length - 1]) <= 0.05 && Integer.parseInt(ss[0]) > 1) {
                for (int i = 2; i <= Integer.parseInt(ss[0]); i++) {
                    path = path + "\t" + ss[i];
                }
                if (max_len2 < Integer.parseInt(ss[0])) {
                    max_len2 = Integer.parseInt(ss[0]);
                }
                f2_paths.put(path, s);
            }
        }

        for (String t : f1_paths.keySet()) {
            ss = t.split("\t");
            tmp = "\t";

            if (f2_paths.containsKey(t)) {
                count++;
                wr.write(f1_paths.get(t) + "\t" + f2_paths.get(t) + "\n");
            }
        }
        System.out.println("size f1 " + f1_paths.size() + " size f2 " + f2_paths.size());
        rd1.close();
        rd2.close();
        wr.close();
    }

    /**
     * Compare two files with nodes
     *
     * @param all Network
     * @param hugo_by_id HGNC ids map
     * @param f1 File 1
     * @param f2 File 2
     * @param out Output file
     */
    public void compare_two_nodes_files(Network all, Map<String, String[]> hugo_by_id, String f1, String f2, String out) throws FileNotFoundException, IOException {
        BufferedReader rd1 = new BufferedReader(new FileReader(f1));
        BufferedReader rd2 = new BufferedReader(new FileReader(f2));
        BufferedWriter wr = new BufferedWriter(new FileWriter(out));
        String s;
        String[] ss;
        Map<String, String> f1_hubs = new HashMap();
        Map<String, String> f2_hubs = new HashMap();
        int count = 0;
        while ((s = rd1.readLine()) != null) {
            ss = s.split("\t");
            if (Integer.parseInt(ss[ss.length - 2]) > 1 && Float.parseFloat(ss[ss.length - 1]) <= 0.05) {
                f1_hubs.put(ss[0], s);
            }
        }
        while ((s = rd2.readLine()) != null) {
            ss = s.split("\t");
            if (Integer.parseInt(ss[ss.length - 2]) > 1 && Float.parseFloat(ss[ss.length - 1]) <= 0.05) {
                f2_hubs.put(ss[0], s);
            }
        }

        for (String t : f1_hubs.keySet()) {
            if (f2_hubs.containsKey(t)) {
                count++;
                wr.write(f1_hubs.get(t) + "\t" + f2_hubs.get(t) + "\n");
            }
        }
        System.out.println("size f1 " + f1_hubs.size() + " size f2 " + f2_hubs.size());
        rd1.close();
        rd2.close();
        wr.close();
    }

    /**
     * Compare a file with paths with a hit list
     *
     * @param all Network
     * @param nutils NetworkManager object
     * @param hugo_by_id HGNC ids map
     * @param f1 File with paths
     * @param out OUtput file
     */
    public void find_hitgenes_on_paths(Network all, NetworkManager nutils, Map<String, String[]> hugo_by_id, String f1, String out) throws FileNotFoundException, IOException {
        BufferedReader rd1 = new BufferedReader(new FileReader(f1));
        BufferedWriter wr = new BufferedWriter(new FileWriter(out));
        String s, path;
        String[] ss, tt;
        Map<String, String> f1_paths = new HashMap();

        int max_len1 = 0, max_len2 = 0;
        int count = 0;
        while ((s = rd1.readLine()) != null) {
            ss = s.split("\t");
            if (Float.parseFloat(ss[ss.length - 1]) <= 0.05 && Integer.parseInt(ss[0]) > 1) {
                path = ss[1];
                for (int i = 2; i <= Integer.parseInt(ss[0]); i++) {
                    path = path + "\t" + ss[i];
                }

                f1_paths.put(path, s);

            }
        }

        for (String w : nutils.hg.keySet()) {
            //System.out.println(w + "\t" + nutils.hg.get(w)+"\n");
        }

        for (String t : f1_paths.keySet()) {
            ss = t.split("\t");
            for (String q : ss) {
                //  System.out.println(all.interactions.get(q).int1.id);
                //  System.out.println(all.interactions.get(q).int2.id);

                if ((nutils.hg.containsKey(all.interactions.get(q).int1.id) || nutils.hg.containsKey(all.interactions.get(q).int2.id)) // && !all.interactions.get(q).int1.id.equals("HGNC:129")
                        // && !all.interactions.get(q).int1.id.equals("pHGNC:129")
                        // && !all.interactions.get(q).int2.id.equals("HGNC:129")
                        // && !all.interactions.get(q).int2.id.equals("pHGNC:129")
                        ) {
                    count++;
                    wr.write(f1_paths.get(t) + "\n");
                    break;
                }
            }

        }
        System.out.println("size f1 " + f1_paths.size());
        rd1.close();
        wr.close();
    }

    private void get_expression(String exp, String gene_list, String out) throws FileNotFoundException, IOException {
        System.out.println("+++++++++++++get_expression++++++++++++");
        BufferedReader rd = new BufferedReader(new FileReader(exp));
        BufferedWriter wr = new BufferedWriter(new FileWriter(out));
        Map<String, List<String>> expr = new HashMap();
        List<String> lst;
        String s;
        String[] ss;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            if (expr.containsKey(ss[1])) {
                expr.get(ss[1]).add(ss[0]);
            } else {
                lst = new ArrayList();
                lst.add(ss[0]);
                expr.put(ss[1], lst);
            }
        }
        rd.close();
        rd = new BufferedReader(new FileReader(gene_list));
        while ((s = rd.readLine()) != null) {
            wr.write(s);
            if (expr.containsKey(s)) {
                for (String t : expr.get(s)) {
                    wr.write("\t" + t);
                }
            }
            wr.newLine();
        }
        expr = null;
        rd.close();
        wr.close();
    }

    /**
     *
     * @param all the value of all
     * @param f_pathways the value of f_pathways
     * @param gw_file the value of gw_file
     * @param conv_table the value of conv_table
     * @param out the value of out
     * @param hugo_by_id the value of hugo_by_id
     * @param centralityManager the value of centralityManager
     * @throws IOException
     */
    public void add_aggregated_pvalues_to_pathways(Network all, String f_pathways, String gw_file, String conv_table, String out, Map<String, String[]> hugo_by_id, CentralityManager centralityManager) throws IOException {
        BufferedReader rd = new BufferedReader(new FileReader(f_pathways));
        BufferedReader rd_pv = new BufferedReader(new FileReader(gw_file));
        BufferedWriter wr = new BufferedWriter(new FileWriter(out));
        int pos = 12;
        String s;
        String[] ss;
        String[] tt;
        int len;
        int number;
        double pvalue;
        double tfct;
        double tfct_test;
        double foldchange;
        double[] tfct_permuted = new double[1000];
        Map<String, List<Float>> permutation_table = new HashMap();
        Map<String, Float> p_values = new HashMap();
        Map<String, Float> fold_changes = new HashMap();
        List<String> genes;
        List<String> genes_new;
        List<String> genes_hugo;
        permutation_table = centralityManager.return_permutation_distribution(gw_file, conv_table, 1000);
        Map<String, Float[]> table = new HashMap();
        table = centralityManager.load_screening_data(gw_file, conv_table);
        for (String q : table.keySet()) {
            p_values.put(q, table.get(q)[1]);
            fold_changes.put(q, table.get(q)[0]);
        }
        String pth;
        String[] pth_ss;
        int pth_len;
        while ((s = rd.readLine()) != null) {
            pth = s;
            pth_ss = pth.split("\t");
            pth_len = pth_ss.length;
            genes_hugo = new ArrayList();
            for (int i = 2; i < pth_len - 1; i++) {
                try {
                    genes_hugo.add(all.interactions.get(pth_ss[i]).int1.id);
                } catch (NullPointerException e) {
                    System.out.println(pth_ss[i]);
                    System.out.println(all.interactions.get(pth_ss[i]).sourcedbentry.get(0));
                    break;
                }
                genes_hugo.add(all.interactions.get(pth_ss[i]).int2.id);
            }
            String symbol = "";
            genes = new ArrayList();
            genes_new = new ArrayList();
            for (String gene : genes_hugo) {
                try {
                    if (gene.startsWith("p")) {
                        symbol = hugo_by_id.get(gene.substring(1))[1];
                    } else {
                        symbol = hugo_by_id.get(gene)[1];
                    }
                } catch (NullPointerException e) {
                    symbol = gene;
                }
                if (!genes.contains(symbol)) {
                    genes.add(symbol);
                }
            }
            tfct = 0.0;
            foldchange = 0.0;
            for (String gene : genes) {
                if (p_values.containsKey(gene)) {
                    genes_new.add(gene);
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
            }
            for (int i = 0; i < 1000; i++) {
                tfct_test = 0.0;
                for (String gene : genes_new) {
                    tfct_test = tfct_test + (Math.log(permutation_table.get(gene).get(i)) * 1.0) / Math.log(2);
                }
                tfct_test = -2 * tfct_test;
                if (tfct_test >= tfct) {
                    number++;
                }
                tfct_permuted[i] = tfct_test;
            }
            pvalue = (number * 1.0) / 1000;
            pth_len = pth_len - 3;
            wr.write(foldchange + "\t" + pvalue + "\t" + pth_len + "\t" + s + "\n");
        }
        wr.close();
        rd.close();
        rd_pv.close();
    }

}

class ListComparator implements Comparator<String[]> {

    @Override
    public int compare(String[] a, String[] b) {
        return Integer.parseInt(a[0]) > Integer.parseInt(b[0]) ? -1 : Integer.parseInt(a[0]) == Integer.parseInt(b[0]) ? 0 : 1;
    }
}

class String_Int {

    String str;
    Integer it;

    public String_Int(String str, Integer it) {
        this.str = str;
        this.it = it;
    }
}

class MaxOccComparator implements Comparator<String_Int> {

    @Override
    public int compare(String_Int a, String_Int b) {
        return a.it > b.it ? -1 : b.it == a.it ? 0 : 1;
    }
}

class EdgeAttributes {

    String beg = "";
    String end = "";
    String[] used_fp = new String[]{"0", "0", "0", "0"};
    Integer overall_occ = 0;
    int[] detailed_occ = new int[]{0, 0, 0, 0};

    public EdgeAttributes(String[] used_fp, int[] detailed_occ, int overall_occ, String beg, String end) {
        this.detailed_occ = detailed_occ;
        this.overall_occ = overall_occ;
        this.used_fp = used_fp;
        this.beg = beg;
        this.end = end;
    }

    public void setAttributes(String[] used_fp, int[] detailed_occ, int overall_occ) {
        this.detailed_occ = detailed_occ;
        this.overall_occ = overall_occ;
        this.used_fp = used_fp;
    }
}

class Link {

    //int path;
    String path_id;
    int offset;
    int length;

    //public Link(int path, int offset, int length) {
    public Link(String path_id, int offset, int length) {
        //this.path = path;
        this.path_id = path_id;
        this.offset = offset;
        this.length = length;
    }
}

class AuxCytoClass {

    String id;
    int numOcc;
    int length;
    List<String> intr;
}

class Path_strength {

    String path;
    Strength strength;

    Path_strength(String path, Strength strength) {
        this.path = path;
        this.strength = strength;
    }
}

class Strength {

    String type;
    int len;

    Strength(String type, int len) {
        this.type = type;
        this.len = len;
    }
}

/*
 public void find_the_shortest_paths(String foundf, String outname, Network all, Map<String, String> hg, Map<String, String> fpl, int d_ppi, int d_tf, int d_mirna, int d_kegg) throws FileNotFoundException, IOException {
 System.out.println("+++++++++++++find_the_shortest_paths++++++++++++");
 BufferedReader rd = new BufferedReader(new FileReader(foundf));
 BufferedWriter wr = new BufferedWriter(new FileWriter(foundf + outname));
 Map<String, List<Path_strength>> hps = new HashMap(); //hg - path- strength
 List<Path_strength> tmp_list;
 Map<String, int[]> min_length = new HashMap();
 int[] tmp_min_len;
 // 0 - ppi; 1- tf; 2- minra; 3 -kegg;
 String s, type, tmp;
 String[] ss;
 while ((s = rd.readLine()) != null) {
 ss = s.split("\t");
 tmp = "";
 if (s.contains("KEGG")) {
 type = "kegg";
 } else {
 if (s.contains("transmir") && s.contains("MIRT")) {
 type = "mirna";
 } else {
 for (int i = 3; i < ss.length; i++) {
 tmp = tmp + "\t" + ss[i];
 }
 if (tmp.contains("man") || tmp.contains("tfacts:") || tmp.contains("XN000")) {
 type = "tf";
 } else {
 type = "ppi";
 }
 }
 }
 if (hps.containsKey(ss[1])) {
 tmp_list = hps.get(ss[1]);
 tmp_list.add(new Path_strength(s, new Strength(type, ss.length - 3)));
 hps.put(ss[1], tmp_list);
 if (type.equals("ppi") && ss.length - 3 < min_length.get(ss[1])[0]) {
 tmp_min_len = min_length.get(ss[1]);
 tmp_min_len[0] = ss.length - 3;
 min_length.put(ss[1], tmp_min_len);
 }
 try {
 if (type.equals("tf") && ss.length - 3 < min_length.get(ss[1])[1]) {
 tmp_min_len = min_length.get(ss[1]);
 tmp_min_len[1] = ss.length - 3;
 min_length.put(ss[1], tmp_min_len);
 }
 } catch (NullPointerException e) {
 System.out.println(s);
 }
 if (type.equals("mirna") && ss.length - 3 < min_length.get(ss[1])[2]) {
 tmp_min_len = min_length.get(ss[1]);
 tmp_min_len[2] = ss.length - 3;
 min_length.put(ss[1], tmp_min_len);
 }
 if (type.equals("kegg") && ss.length - 3 < min_length.get(ss[1])[3]) {
 tmp_min_len = min_length.get(ss[1]);
 tmp_min_len[3] = ss.length - 3;
 min_length.put(ss[1], tmp_min_len);
 }
 } else {
 tmp_list = new ArrayList();
 tmp_list.add(new Path_strength(s, new Strength(type, ss.length - 3)));
 hps.put(ss[1], tmp_list);

 tmp_min_len = new int[]{Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE};
 if (type.equals("ppi")) {
 tmp_min_len[0] = ss.length - 3;
 min_length.put(ss[1], tmp_min_len);
 }
 if (type.equals("tf")) {
 tmp_min_len[1] = ss.length - 3;
 min_length.put(ss[1], tmp_min_len);
 }
 if (type.equals("mirna")) {
 tmp_min_len[2] = ss.length - 3;
 min_length.put(ss[1], tmp_min_len);
 }
 if (type.equals("kegg")) {
 tmp_min_len[3] = ss.length - 3;
 min_length.put(ss[1], tmp_min_len);
 }
 }
 }
 for (String t : hps.keySet()) {
 System.out.println(t + "\t" + min_length.get(t)[0] + "\t" + min_length.get(t)[1] + "\t" + min_length.get(t)[2] + "\t" + min_length.get(t)[3]);
 for (Path_strength ps : hps.get(t)) {
 if (ps.strength.type.equals("ppi") & ps.strength.len <= min_length.get(t)[0] + d_ppi) {
 wr.write(ps.path + "\n");
 }
 if (ps.strength.type.equals("tf") & ps.strength.len <= min_length.get(t)[1] + d_tf) {
 wr.write(ps.path + "\n");
 }
 if (ps.strength.type.equals("mirna") & ps.strength.len <= min_length.get(t)[2] + d_mirna) {
 wr.write(ps.path + "\n");
 }
 if (ps.strength.type.equals("kegg") & ps.strength.len <= min_length.get(t)[3] + d_kegg) {
 wr.write(ps.path + "\n");
 }
 }
 }
 rd.close();
 wr.close();
 }


 public void getPathways(String foundf, Network all, Map<String, String> hg, Map<String, String> fpl) throws FileNotFoundException, IOException {
 BufferedReader rd = new BufferedReader(new FileReader(foundf));
 BufferedReader rd2 = new BufferedReader(new FileReader("F:\\Dropbox\\_Работа\\Programms\\DATA\\Muscle Differentiation\\miRNA results\\JAG1_all_paths.txt"));
 String s;
 String[] ss, ss_tmp;
 String hit = "pHGNC:6188";
 for (String t : hg.keySet()) {
 //System.out.println(t + " " + hg.get(t));
 }
 while ((s = rd.readLine()) != null) {
 ss = s.split("\t");
 ss_tmp = ss[ss.length - 1].split(" ");
 if (ss_tmp[ss_tmp.length - 1].equals(hit) && s.contains("MIRT")) {
 // System.out.println(s); // results in file 
 }
 }
 int count = 0;
 int id_path = 0;
 String[] fpls = {"pHGNC:7611", "pHGNC:4223", "pHGNC:7566", "pHGNC:5466"};
 Map<String, Integer> fpl_count = new HashMap();
 fpl_count.put("pHGNC:7611", 0);
 fpl_count.put("pHGNC:4223", 0);
 fpl_count.put("pHGNC:7566", 0);
 fpl_count.put("pHGNC:5466", 0);
 while ((s = rd2.readLine()) != null) {
 ss = s.split("\t");
 count = fpl_count.get(ss[0]);
 fpl_count.put(ss[0], count + 1);
 for (int i = 1; i < ss.length - 1; i++) {
 System.out.println(ss[0] + "_" + fpl_count.get(ss[0]) + "_" + i + "\t" + all.interactions.get(ss[i]).int1.id + "\t" + all.interactions.get(ss[i]).int2.id + "\t" + ss[0] + "\t" + ss[0] + "_" + fpl_count.get(ss[0]) + "\t" + ss[i]);
 }
 }
 rd.close();
 rd2.close();
 }


 /*  List<Map<String, Integer>> arg_l;
 Map<String, Integer> arg_ll;
 if (mi_list.containsKey(mirtarbase.interactions.get(t).int1.id)) {
 arg_l = mi_list.get(mirtarbase.interactions.get(t).int1.id);
 for (int i = 0; i < 3; i++) {
 arg_ll = arg_l.get(i);
 for (String tt : id_list.get(t).get(i).keySet()) {
 if (arg_ll.containsKey(tt)) {
 d = arg_ll.get(tt) + id_list.get(t).get(i).get(tt);
 ln_m_tmp.put(tt, d);
 } else {
 arg_ll.put(tt, id_list.get(t).get(i).get(tt));
 }
 }
 arg_l.set(i, arg_ll);
 }

 mi_list.put(mirtarbase.interactions.get(t).int1.id, arg_l);
 } else {
 mi_list.put(mirtarbase.interactions.get(t).int1.id, id_list.get(t));
 } 






 *1
 beg = fpl.get(ss_tmp[0]);
 // ss_tmp = ss[ss.length - 1].split(" ");
 //end = fpl.get(ss_tmp[ss_tmp.length - 1]);                    
 //end = hg.get(ss_tmp[ss_tmp.length - 1]);
 //ss_tmp = ss[0].split(" ");
 //beg = hg.get(ss_tmp[0]);

 *2
 beg = fpl.get(ss_tmp[0]);
 // ss_tmp = ss[ss.length - 1].split(" ");
 //end = fpl.get(ss_tmp[ss_tmp.length - 1]);                    
 //end = hg.get(ss_tmp[ss_tmp.length - 1]);
 //ss_tmp = ss[0].split(" ");
 //beg = hg.get(ss_tmp[0]);
 */
/*




 public void create_network_for_cytoscape(String foundf, Network all, Map<String, String> hg, Map<String, String> fpl) throws FileNotFoundException, IOException {
 BufferedReader rd = new BufferedReader(new FileReader(foundf));
 BufferedWriter wr = new BufferedWriter(new FileWriter(foundf + "_cyto"));
 String s,fp;
 String[] ss;
 int count = 0;
 int id_path = 0;
 int ss_len;
 String[] fpls = {"pHGNC:7611", "pHGNC:4223", "pHGNC:7566", "pHGNC:5466"};
 Map<String, Integer> fpl_count = new HashMap();
 fpl_count.put("pHGNC:7611", 0);
 fpl_count.put("pHGNC:4223", 0);
 fpl_count.put("pHGNC:7566", 0);
 fpl_count.put("pHGNC:5466", 0);
        
 while ((s = rd.readLine()) != null) {
 ss = s.split("\t");
 ss_len = ss.length;
 fp = ss[ss_len - 1];
 count = fpl_count.get(fp);
 fpl_count.put(ss[0], count + 1);
            
 for (int i = 1; i < ss.length - 1; i++) {
 wr.write(fp + "_" + fpl_count.get(fp) + "_" + i + "\t" + all.interactions.get(ss[i]).int1.id + "\t" + all.interactions.get(ss[i]).int2.id + "\t" + ss[ss_len - 1] + "\t" + ss[ss_len - 1] + "_" + fpl_count.get(fp) + "\t" + ss[i] + "\n");
 }
 }
 rd.close();
 wr.close();
 }



 public void getPathways2_JAG1_to_delete() throws FileNotFoundException, IOException {
 BufferedReader rd = new BufferedReader(new FileReader(path + jag1_in));
 //BufferedWriter wr = new BufferedWriter(new FileWriter("F:\\Dropbox\\_Работа\\Programms\\DATA\\Muscle Differentiation\\miRNA results\\JAG1_all_paths_interactions_edt.txt"));
 BufferedWriter wr2 = new BufferedWriter(new FileWriter(path + jag1_out));
 String s, beg, end;
 String[] ss, tmp;
 int[] arr;
 int count;
 //Map<String, String[]> freq = new HashMap();
 Map<String, EdgeAttributes> occur = new HashMap();
 // pHGNC:7611 , pHGNC:4223 , pHGNC:7566  , pHGNC:5466

 String[] fp = new String[]{"pHGNC:7611", "pHGNC:4223", "pHGNC:7566", "pHGNC:5466"};
 while ((s = rd.readLine()) != null) {
 ss = s.split("\t");

 if (occur.containsKey(ss[5])) {
 tmp = occur.get(ss[5]).used_fp;
 arr = occur.get(ss[5]).detailed_occ;
 count = occur.get(ss[5]).overall_occ;
 for (int i = 0; i < fp.length; i++) {
 if (ss[3].equals(fp[i])) {
 tmp[i] = "1";
 arr[i]++;
 count++;
 }
 }
 occur.get(ss[5]).setAttributes(tmp, arr, count);
 occur.put(ss[5], occur.get(ss[5]));
 } else {
 tmp = new String[]{"0", "0", "0", "0"};
 arr = new int[]{0, 0, 0, 0};
 count = 0;
 beg = "";
 end = "";
 for (int i = 0; i < fp.length; i++) {
 if (ss[3].equals(fp[i])) {
 tmp[i] = "1";
 arr[i]++;
 count++;
 beg = ss[1];
 end = ss[2];
 }
 }
 occur.put(ss[5], new EdgeAttributes(tmp, arr, count, beg, end));
 }
 // wr.write(s + "\t" + "OCC" + tmp[0] + tmp[1] + tmp[2] + tmp[3] );
 // wr.newLine();            
 }
 rd.close();
 rd = new BufferedReader(new FileReader(path + jag1_in));
 for (String t : occur.keySet()) {
 wr2.write(occur.get(t).beg + " (pp) " + occur.get(t).end);
 wr2.write("\t" + "OCC" + occur.get(t).used_fp[0] + occur.get(t).used_fp[1] + occur.get(t).used_fp[2] + occur.get(t).used_fp[3]);
 wr2.write("\t" + occur.get(t).overall_occ);
 wr2.write("\t" + occur.get(t).detailed_occ[0] + "," + occur.get(t).detailed_occ[1] + "," + occur.get(t).detailed_occ[2] + "," + occur.get(t).detailed_occ[3]);
 wr2.newLine();
 }
 rd.close();
 // wr.close();
 wr2.close();
 }

 public void calc_to_delete() throws IOException {
 BufferedReader rd = new BufferedReader(new FileReader(path + fname));
 BufferedWriter wr = new BufferedWriter(new FileWriter(path + outfile));
 String s;
 int i;
 while ((s = rd.readLine()) != null) {
 if (!inter.containsKey(s)) {
 inter.put(s, new Integer(1));
 } else {
 i = 1 + inter.get(s).intValue();
 inter.put(s, new Integer(i));
 }
 }

 for (String ss : inter.keySet()) {
 wr.write(ss + "\t" + inter.get(ss));
 wr.newLine();
 }
 rd.close();
 wr.close();
 }


 private void checkLinks(Map<Integer, String[]> paths, Map<Integer, List<Link>> links, String[] tmp, int path, int offset, int length) {
 boolean eq;
 List<Link> lk, list;
 boolean analyzed = false;
 for (Integer it : links.keySet()) {
 if (analyzed) {
 break;
 }
 eq = true;
 lk = links.get(it);
 if (lk.get(0).length == length) {
 // System.out.println("bbb");
 for (int i = 0; i < length; i++) {
 //System.out.println("i = " + i + " i-offset = " + (i - offset) + " path number " + it + " path size " + paths.get(lk.get(0).path).length);
 if (!paths.get(lk.get(0).path)[lk.get(0).offset + i].equals(tmp[i])) {
 eq = false;
 }
 }
 if (eq) {
 //  System.out.println("ccc");
 lk.add(new Link(path, offset, length));
 links.put(it, lk);
 analyzed = true;
 break;
 }
 }
 }
 if (!analyzed) {
 // System.out.println("dd");
 list = new ArrayList();
 list.add(new Link(path, offset, length));
 links.put(links.size() + 1, list);
 }
 }

 private void writeLinks(Map<Integer, String[]> paths, Map<Integer, List<Link>> links, String outf) throws IOException {
 List<Link> lk;
 BufferedWriter wr = new BufferedWriter(new FileWriter(outf));
 String tmp;
 for (Integer it : links.keySet()) {
 lk = links.get(it);
 tmp = lk.size() + "\t" + lk.get(0).length + "\t" + paths.get(lk.get(0).path)[0] + "\t" + paths.get(lk.get(0).path)[paths.get(lk.get(0).path).length - 1];
 for (int i = lk.get(0).offset; i < lk.get(0).length + lk.get(0).offset; i++) {
 tmp = tmp + "\t" + paths.get(lk.get(0).path)[i];
 }
 wr.write(tmp);
 wr.newLine();
 }
 wr.close();
 }

 public void find_overreprs(int max_length, int min_length, String inf, String outf) throws IOException {
 //System.out.println("+++++++++++++analysePaths++++++++++++");
 BufferedReader rd = new BufferedReader(new FileReader(inf));
 Map<Integer, List<Link>> links = new HashMap();
 Map<Integer, String[]> paths = new HashMap();
 String s;
 String[] ss, tmpa;
 //String tmp = "";
 int min, max, max_end;//, tmpi;
 int ii = 0;
 int u = 0;
 while ((s = rd.readLine()) != null) {
 ss = s.split("\t");
 paths.put(ii, ss);
 ii++;
 }
 int count = 0;
 for (Integer it : paths.keySet()) {
 count++;
 if (count >= 10000) {
 break;
 }
 ss = paths.get(it);
 if (ss.length - 2 < max_length) {
 max = ss.length - 2;
 } else {
 max = max_length;
 }
 min = min_length;
 for (int i = 1; i <= ss.length - 1 - max; i++) {
 for (int j = min; j <= max; j++) {
 tmpa = new String[j];
 u = 0;
 for (int k = i; k < i + j; k++) {
 tmpa[u] = ss[k];
 //  tmp = tmp + "\t" + ss[k];
 u++;
 }
 checkLinks(paths, links, tmpa, it, i, j);
 }
 }
 max_end = max;
 for (int i = ss.length - 1 - max + 1; i <= ss.length - 1 - min; i++) {
 max_end = max_end - 1;
 for (int j = min; j <= max_end; j++) {
 tmpa = new String[j];
 u = 0;
 for (int k = i; k < i + j; k++) {
 tmpa[u] = ss[k];
 // tmp = tmp + "\t" + ss[k];
 u++;
 }
 checkLinks(paths, links, tmpa, it, i, j);
 }
 }
 }
 writeLinks(paths, links, outf);
 //System.out.println("links size " + links.size());
 }


 public void find_miRNAs_on_pathways(String foundf, String resf, Network mirtarbase, Map<String, String> hg, Map<String, String> fpl, int length, int min, String mask) throws FileNotFoundException, IOException {
 System.out.println("+++++++++++++find_miRNAs_on_pathways++++++++++++");
 BufferedReader rd = new BufferedReader(new FileReader(foundf));
 //BufferedWriter wr= new BufferedWriter (new FileWriter(resf));
 String s;//, mi;
 String[] ss;
 int count;
 Map<String, Integer> id_count = new HashMap();
 Map<String, int[]> id_count_by_3 = new HashMap();        //mirtarbase id - data
 Map<String, List<Map<String, Integer>>> id_list = new HashMap();

 Map<String, Integer> mi_count = new HashMap(); //mirna - number of occurrence
 Map<String, int[]> mi_count_by_3 = new HashMap();        //mirna - data
 Map<String, List<Map<String, Integer>>> mi_list = new HashMap();

 //int[] ln = new int[3];
 int[] ln = new int[length - min];
 List<Map<String, Integer>> ln_m = new ArrayList();
 Map<String, Integer> ln_m_tmp = new HashMap();
 String beg, end;
 String[] ss_tmp;
 List<String> all_mi;
 int d;
 int start_pos;
 while ((s = rd.readLine()) != null) {
 if (s.contains("MIRT")) {
 all_mi = new ArrayList();
 ss = s.split("\t");
 if (mask.equals("countonlymiddle")) {
 start_pos = 2;
 } else {
 start_pos = 0;
 }
 for (int i = start_pos; i < ss.length; i++) {
 if (ss[i].contains("MIRT")) {
 all_mi.add(ss[i]);
 }
 }

 //  mi = s.substring(s.indexOf("MIRT"), s.indexOf("MIRT") + 10);
 for (String mi : all_mi) {
 if (id_count.containsKey(mi)) {
 ////////////////////////////////
 count = id_count.get(mi);
 id_count.put(mi, count + 1);
 ////////////////////////////////
 ss = s.split("\t");
 ln = id_count_by_3.get(mi);
 for (int i = 0; i < length - min; i++) {
 if (ss.length == i + 5) {
 ln[i]++;
 }
 }
 id_count_by_3.put(mi, ln);

 ////////////////////////////////
 //*1
 beg = hg.get(ss[0]);
 end = fpl.get(ss[ss.length - 1]);
 //System.out.println("-----------------------------------------------beg-end " +mi + " " + beg + "-" + end);
 ln_m = id_list.get(mi);
 for (int i = 0; i < length - min; i++) {
 if (ss.length == i + 5) {
 ln_m_tmp = ln_m.get(i);
 if (ln_m.get(i).containsKey(beg + "-" + end)) {
 d = ln_m_tmp.get(beg + "-" + end) + 1;
 ln_m_tmp.put(beg + "-" + end, d);
 } else {
 ln_m_tmp.put(beg + "-" + end, 1);
 }
 ln_m.set(i, ln_m_tmp);
 }
 }
 id_list.put(mi, ln_m);
 ////////////////////////////////
 } else {
 ////////////////////////////////
 id_count.put(mi, 1);
 //////////////////////////////
 ln = new int[length - min];
 for (int i = 0; i < length - min; i++) {
 ln[i] = 0;
 }
 ss = s.split("\t");
 for (int i = 0; i < length - min; i++) {
 ln[i] = 0;
 if (ss.length == i + 5) {
 //  mirna_int_ln.get(mi)[i]++;
 ln[i] = 1;
 }
 }
 id_count_by_3.put(mi, ln);
 ////////////////////////////////
 //*2
 beg = hg.get(ss[0]);  // bez + beg\end
 end = fpl.get(ss[ss.length - 1]);
 //System.out.println("-----------------------------------------------beg-end" + beg + "-" + end);
 ln_m = new ArrayList();
 for (int i = 0; i < length - min; i++) {
 ln_m.add(new HashMap());
 }
 ln_m_tmp = new HashMap();
 for (int i = 0; i < length - min; i++) {
 if (ss.length == i + 5) {
 ln_m_tmp.put(beg + "-" + end, 1);
 ln_m.set(i, ln_m_tmp);
 }
 }
 id_list.put(mi, ln_m);
 }
 }
 }
 }
 count = 0;
 int cc = 0;
 int[] arg;

 int check1 = 0;
 int check2 = 0;
 for (String t : id_count.keySet()) {
 count = count + id_count.get(t);
 if (mi_count.containsKey(mirtarbase.interactions.get(t).int1.id)) {
 cc = mi_count.get(mirtarbase.interactions.get(t).int1.id);
 mi_count.put(mirtarbase.interactions.get(t).int1.id, cc + id_count.get(t));
 } else {
 mi_count.put(mirtarbase.interactions.get(t).int1.id, id_count.get(t));
 }

 ///////////////////////////////////////////////////////////////////////////////////
 if (mi_count_by_3.containsKey(mirtarbase.interactions.get(t).int1.id)) {
 arg = mi_count_by_3.get(mirtarbase.interactions.get(t).int1.id);
 for (int i = 0; i < length - min; i++) {
 arg[i] = arg[i] + id_count_by_3.get(t)[i];
 }
                
 mi_count_by_3.put(mirtarbase.interactions.get(t).int1.id, arg);
 } else {
 mi_count_by_3.put(mirtarbase.interactions.get(t).int1.id, id_count_by_3.get(t));
 }

 ///////////////////////////////////////////////////////////////////////////////////
 //**
 ///////////////////////////////////////////////////////////////////////////////////
 if (mirtarbase.interactions.get(t).int1.id.equals("hsa-mir-1")) {
 check1 = check1 + id_count.get(t);
 for (int i = 0; i < length - min; i++) {
 for (String tt : id_list.get(t).get(i).keySet()) {
 check2 = check2 + id_list.get(t).get(i).get(tt);
 System.out.println("--- " + tt + " " + id_list.get(t).get(i).get(tt));
 }
 }
 }
 ///////////////////////////////////////////////////////////////////////////////////
 }
 String mir = "";
 //Map<String, List<Map<String, Integer>>> id_list = new HashMap();
 List<Map<String, Integer>> tmp_list = new ArrayList();
 Map<String, Integer> tmp_map = new HashMap();
 for (String t : id_list.keySet()) {
 d = 0;
 mir = mirtarbase.interactions.get(t).int1.id;
 if (mi_list.containsKey(mir)) {
 for (int i = 0; i < length - min; i++) {
 for (String tt : id_list.get(t).get(i).keySet()) {
 tmp_list = new ArrayList();
 tmp_map = new HashMap();
 if (mi_list.get(mir).get(i).containsKey(tt)) {
 d = mi_list.get(mir).get(i).get(tt) + id_list.get(t).get(i).get(tt);
 tmp_list.addAll(mi_list.get(mir));
 //tmp_map=id_list.get(t).get(i);
 tmp_map.putAll(id_list.get(t).get(i));
 tmp_map.put(tt, d);
 tmp_list.set(i, tmp_map);
 // mi_list.put(mir, tmp_list);
 mi_list.get(mir).get(i).put(tt, d);
 } else {
 d = id_list.get(t).get(i).get(tt);
 tmp_list.addAll(mi_list.get(mir));
 //tmp_map=id_list.get(t).get(i);
 tmp_map.putAll(id_list.get(t).get(i));
 tmp_map.put(tt, d);
 tmp_list.set(i, tmp_map);
 //  mi_list.put(mir, tmp_list); 
 mi_list.get(mir).get(i).put(tt, id_list.get(t).get(i).get(tt));
 }
 //mi_list.get(mirtarbase.interactions.get(t).int1.id).get(i) 
 }
 }
 } else {
 mi_list.put(mir, id_list.get(t));
 }
 }
 System.out.println("check1 " + check1);
 System.out.println("check2 " + check2);
 System.out.println("count " + count);
 count = 0;
 for (String t : mi_count.keySet()) {
 count = count + mi_count.get(t);
 System.out.print(t + "\t" + mi_count.get(t) + "\t");
 for (int i = 0; i < length - min; i++) {
 System.out.print(mi_count_by_3.get(t)[i] + "\t");
 }
 //System.out.print(t + "\t" + mi_count.get(t) + "\t" + mi_count_by_3.get(t)[0] + "\t" + mi_count_by_3.get(t)[1] + "\t" + mi_count_by_3.get(t)[2] + "\t");
 for (int i = 0; i < length - min; i++) {
 //  System.out.print(i + 3 + ":\t");
 System.out.print(i + length - min - 1 + ":\t");
 for (String tt : mi_list.get(t).get(i).keySet()) {
 System.out.print(tt + " x " + mi_list.get(t).get(i).get(tt) + "\t");
 }
 }
 System.out.println();
 }
 System.out.println("count " + count);
 rd.close();
 }

 public void create_network_for_cytoscape_only_best(String foundf_cyto, Network all, Map<String, String> hg, Map<String, String> fpl) throws FileNotFoundException, IOException {
 BufferedReader rd = new BufferedReader(new FileReader(foundf_cyto));
 BufferedWriter wr = new BufferedWriter(new FileWriter(foundf_cyto + "_cyto_only_best"));
 String s, fp;
 String[] ss;
 int count = 0;
 int id_path = 0;
 int ss_len;
 Map<String, AuxCytoClass> intr = new HashMap();
 boolean checked = false;
 while ((s = rd.readLine()) != null) {
 ss = s.split("\t");
 checked = false;
 ss = s.split("\t");
 for (String t : intr.keySet()) {
 if (intr.get(t).id.equals(ss[5])) {
 if (intr.get(t).length > Integer.getInteger(ss[4]) && intr.get(t).numOcc > Integer.getInteger(ss[3])) {
 //change
 }
 if (intr.get(t).length < Integer.getInteger(ss[4]) && intr.get(t).numOcc < Integer.getInteger(ss[3])) {
 //add
 }
 if ((intr.get(t).length >= Integer.getInteger(ss[4]) && intr.get(t).numOcc <= Integer.getInteger(ss[3])) || (intr.get(t).length <= Integer.getInteger(ss[4]) && intr.get(t).numOcc >= Integer.getInteger(ss[3]))) {
 //keep both  
 }
 }
 }
 /*     for (int i = 4; i < ss.length; i++) {
 wr.write(fp + "_" + fpl_count.get(fp) + "_" + i+ "\t"); 
 wr.write(ss[2]+ "\t");
 wr.write(ss[3]+ "\t");
 wr.write(ss[0]+ "\t");
 wr.write(ss[1]+ "\t");
 wr.write(ss[i]+ "\t");                
 wr.write(all.interactions.get(ss[i]).int1.id+ "\t");
 wr.write(all.interactions.get(ss[i]).int2.id+ "\t");
 wr.write(fp + "_" + fpl_count.get(fp)+ "\t");
 wr.newLine();
 } */
/*     }
 rd.close();
 wr.close();
 }


 public void compare_two_top_ranked(Network all, Map<String, String[]> hugo_by_id, String f1, String f2, String out) throws FileNotFoundException, IOException {
 BufferedReader rd1 = new BufferedReader(new FileReader(f1));
 BufferedReader rd2 = new BufferedReader(new FileReader(f2));
 BufferedWriter wr = new BufferedWriter(new FileWriter(out));
 String s, path, tmp, n1, n2;
 String[] ss, tt;
 Map<String, String> f1_path_rank = new HashMap();
 Map<String, String> f2_path_rank = new HashMap();
 Map<String, String> f1_path_hgfp = new HashMap();
 Map<String, String> f2_path_hgfp = new HashMap();
 int max_len = 0;
 while ((s = rd1.readLine()) != null) {
 ss = s.split("\t");
 path = ss[4];
 for (int i = 5; i < ss.length; i++) {
 path = path + "\t" + ss[i];
 }
 f1_path_rank.put(path, ss[0] + "\t" + ss[1]);
 f1_path_hgfp.put(path, ss[2] + "\t" + ss[3]);
 if (max_len < ss.length - 4) {
 max_len = ss.length - 4;
 }
 }
 while ((s = rd2.readLine()) != null) {
 ss = s.split("\t");
 path = ss[4];
 for (int i = 5; i < ss.length; i++) {
 path = path + "\t" + ss[i];
 }
 f2_path_rank.put(path, ss[0] + "\t" + ss[1]);
 f2_path_hgfp.put(path, ss[2] + "\t" + ss[3]);
 }
 for (String t : f1_path_rank.keySet()) {
 ss = t.split("\t");
 tmp = "\t";
 if (ss.length < max_len) {
 for (int i = 0; i < max_len - ss.length; i++) {
 tmp = tmp + "\t";
 }
 }
 if (f2_path_rank.containsKey(t)) {
 wr.write(f1_path_rank.get(t) + "\t" + f1_path_hgfp.get(t) + "\t" + t + "\t" + tmp + "\t" + f2_path_rank.get(t) + "\t" + f2_path_hgfp.get(t) + "\t");
 tt = t.split("\t");
 for (String tt1 : tt) {
 try {
 if (all.interactions.get(tt1).int1.id.startsWith("hsa-")) {
 n1 = all.interactions.get(tt1).int1.id;
 } else if (all.interactions.get(tt1).int1.id.startsWith("p")) {
 n1 = "p" + hugo_by_id.get(all.interactions.get(tt1).int1.id.substring(1))[1];
 } else {
 n1 = hugo_by_id.get(all.interactions.get(tt1).int1.id)[1];
 }
 } catch (NullPointerException e) {
 n1 = all.interactions.get(tt1).int1.id;
 }
 try {
 if (all.interactions.get(tt1).int2.id.startsWith("hsa-")) {
 n2 = all.interactions.get(tt1).int2.id;
 } else if (all.interactions.get(tt1).int2.id.startsWith("p")) {
 n2 = "p" + hugo_by_id.get(all.interactions.get(tt1).int2.id.substring(1))[1];
 } else {
 n2 = hugo_by_id.get(all.interactions.get(tt1).int2.id)[1];
 }
 } catch (NullPointerException e) {
 n2 = all.interactions.get(tt1).int2.id;
 }
 wr.write(" \t" + n1 + "-" + n2);
 }
 wr.newLine();
 }
 }
 rd1.close();
 rd2.close();
 wr.close();
 }
 */

/*


 String path = "F:\\Dropbox\\_Работа\\Programms\\DATA\\Muscle Differentiation\\";
 String fname = "nw2.txt";
 String outfile = "nw2_edt";

 String jag1_in = "miRNA results\\JAG1_all_paths_interactions.txt";
 String jag1_out = "\\miRNA results\\JAG1_all_paths_interactions_edge_attr.txt";
 Map<String, Integer> inter = new HashMap();

 private void compare_two_lists(String l1, String l2, String out) throws FileNotFoundException, IOException {
 System.out.println("+++++++++++++compare two list++++++++++++");
 BufferedReader rd = new BufferedReader(new FileReader(l1));
 BufferedWriter wr = new BufferedWriter(new FileWriter(out));
 Set<String> lst = new HashSet();
 String s;
 String[] ss;
 while ((s = rd.readLine()) != null) {
 lst.add(s);
 }
 rd.close();
 rd = new BufferedReader(new FileReader(l2));
 while ((s = rd.readLine()) != null) {
 if (lst.contains(s)) {
 wr.write(s);
 wr.newLine();
 }
 }
 lst = null;
 rd.close();
 wr.close();

 }

    

 public void compare_two_top_ranked(Network all, Map<String, String[]> hugo_by_id, String f1, String f2, String out, String prefix) throws FileNotFoundException, IOException {
 BufferedReader rd1 = new BufferedReader(new FileReader(f1));
 BufferedReader rd2 = new BufferedReader(new FileReader(f2));
 BufferedWriter wr = new BufferedWriter(new FileWriter(out));
 String s, path, tmp, n1, n2, names;
 String[] ss, tt;
 Map<String, String> f1_path_rank = new HashMap();
 Map<String, String> f2_path_rank = new HashMap();
 Map<String, String> f1_path_hgfp = new HashMap();
 Map<String, String> f2_path_hgfp = new HashMap();
 Map<String, String> f1_path_names = new HashMap();
 Map<String, String> f2_path_names = new HashMap();
 int max_len = 0;
 int count = 0;
 while ((s = rd1.readLine()) != null) {
 ss = s.split("\t");
 path = ss[5];
 for (int i = 6; i < ss.length - Integer.parseInt(ss[1]); i++) {
 path = path + "\t" + ss[i];
 }
 names = ss[ss.length - Integer.parseInt(ss[1])];
 for (int i = ss.length - Integer.parseInt(ss[1]) + 1; i < ss.length; i++) {
 names = names + "\t" + ss[i];
 }
 f1_path_names.put(path, names);
 f1_path_rank.put(path, ss[0] + "\t" + ss[1]);
 f1_path_hgfp.put(path, ss[2] + "\t" + ss[3] + "\t" + ss[4]);
 if (max_len < Integer.parseInt(ss[1])) {
 max_len = Integer.parseInt(ss[1]);
 }
 }
        
 while ((s = rd2.readLine()) != null) {
 ss = s.split("\t");
 path = ss[5];
 for (int i = 6; i < ss.length - Integer.parseInt(ss[1]); i++) {
 path = path + "\t" + ss[i];
 }
 names = ss[ss.length - Integer.parseInt(ss[1])];
 for (int i = ss.length - Integer.parseInt(ss[1]) + 1; i < ss.length; i++) {
 names = names + "\t" + ss[i];
 }
 f2_path_names.put(path, names);
 f2_path_rank.put(path, ss[0] + "\t" + ss[1]);
 f2_path_hgfp.put(path, ss[2] + "\t" + ss[3] + "\t" + ss[4]);
 }
 for (String t : f1_path_rank.keySet()) {
 ss = t.split("\t");
 tmp = "\t";
 if (ss.length < max_len) {
 for (int i = 0; i < max_len - ss.length; i++) {
 tmp = tmp + "\t";
 }
 }

 if (f2_path_rank.containsKey(t)) {
 count++;
 wr.write(prefix + "_" + count + "\t" + f1_path_rank.get(t) + "\t" + f1_path_hgfp.get(t) + "\t" + t + "\t" + tmp + "\t" + f2_path_rank.get(t) + "\t" + f2_path_hgfp.get(t));
 wr.write("\t" + f1_path_names.get(t));
 wr.newLine();
 }
 }
 rd1.close();
 rd2.close();
 wr.close();
 }

 */
