package masterPATH;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * DBManager class contains methods to load databases
 *
 * @author Natalia Rubanova
 */
public class DBManager {

    Network hprd = new Network();
    Network hprd_invivo = new Network();
    Network pc = new Network();
    Network transpath = new Network();
    Network kegg = new Network();
    Network transmir = new Network();
    Network mirtarbase = new Network();
    Network tfacts = new Network();
    Network hippie = new Network();
    Network hippie_invivo = new Network();
    Network hippie_high = new Network();
    Network hippie_medium_meq = new Network();
    Network hippie_medium_m = new Network();
    Network hippie_medium_m0 = new Network();
    Network hippie_all = new Network();
    Network signor = new Network();
    Network signalink = new Network();

    Map<String, String[]> hugo = new HashMap(); //name - строка
    Map<String, String[]> hugo_by_id = new HashMap(); // id -строка

    Map<String, String> tfacts_to_hugo = new HashMap();
    Map<String, String[]> hugo_entrez = new HashMap();
    Map<String, String[]> entrez_wo_hugo = new HashMap();

    Map<String, String> hprd_to_hugo = new HashMap();
    Map<String, String> transpath_to_hugo = new HashMap();
    Map<String, String> bioDBhprd_to_hugo = new HashMap();
    List<String> lostgenes = new ArrayList();

    /**
     *
     * Constructor
     */
    public DBManager(String hugo_path, String hugo_entrez_path, String entrz_wo_hugo, String hprd_bioDB, String tfacts_crosstable) throws IOException {
        loadHUGO(hugo_path, hugo_entrez_path, entrz_wo_hugo, hprd_bioDB, tfacts_crosstable);
    }

    /**
     * Loads HGNC nomenclature
     */
    public void loadHUGO(String hugo_path, String hugo_entrez_path, String entrz_wo_hugo, String hprd_bioDB, String tfacts_crosstable) throws FileNotFoundException, IOException {
        String s;
        String[] ss;
        BufferedReader rd;
        int i;
        //for (i = 0; i < FoldersPaths.hugo_path.length; i++) {
        for (i = 0; i < 1; i++) {
            rd = new BufferedReader(new FileReader(hugo_path));

            s = rd.readLine();
            while ((s = rd.readLine()) != null) {
                if (!s.toLowerCase().contains("entry withdrawn")) {
                    ss = s.split("\t");
                    //System.out.println(Arrays.toString(ss));
                    hugo.put(ss[1], ss);
                    hugo_by_id.put(ss[0], ss);
                }
            }
            rd.close();
        }

        rd = new BufferedReader(new FileReader(hugo_entrez_path));
        ss = new String[10];
        s = rd.readLine();
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            try {
                hugo_entrez.put(ss[9], ss);
            } catch (ArrayIndexOutOfBoundsException e) {
                // System.out.println(s);
                continue;
            }
        }
        rd.close();

        rd = new BufferedReader(new FileReader(entrz_wo_hugo));
        ss = new String[3];
        s = rd.readLine();
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            entrez_wo_hugo.put(ss[0], ss);
        }
        rd.close();

        rd = new BufferedReader(new FileReader(hprd_bioDB));
        ss = new String[5];
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            bioDBhprd_to_hugo.put(ss[0], ss[4]);
        }
        rd.close();

        rd = new BufferedReader(new FileReader(tfacts_crosstable));
        ss = new String[2];
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            tfacts_to_hugo.put(ss[0], ss[1]);
        }
        rd.close();
    }

    /**
     * Load HPRD database
     */
    public void loadHPRD(String dbpath_conv, String dbpath_int) throws FileNotFoundException, IOException {
        //entities FoldersPaths.hprd_conv FoldersPaths.hprd_int
        System.out.println("+++++++++++++HPRD++++++++++++");
        List<String> db_entr;
        List<String> other_ids = new ArrayList(); //interactions , сначала пустой
        other_ids.add("-");
        String s;
        String[] ss;
        List<String[]> new_ids;
        String new_id;
        BufferedReader rd = new BufferedReader(new FileReader(dbpath_conv));
        HashMap<String, String> conv_map = new HashMap();
        Node n1, n2, n1_invivo, n2_invivo;
        boolean invivo = false;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            conv_map.put(ss[0], ss[ss.length - 1]);
        }
        rd.close();

        //interactions
        rd = new BufferedReader(new FileReader(dbpath_int));
        while ((s = rd.readLine()) != null) {
            if (s.contains("invivo")) {
                invivo = true;
            }
            ss = s.split("\t");
            if (!conv_map.containsKey(ss[2]) || !conv_map.containsKey(ss[5])) {
                continue;
            }
            new_id = "p" + conv_map.get(ss[2]);
            if (!hprd.nodes.containsKey(new_id)) {
                new_ids = new ArrayList();
                new_ids.add(new String[]{ss[2]});
                hprd.nodes.put(new_id, new Node(new_id, "protein", "id", "hprd", new_ids));
                n1 = hprd.nodes.get(new_id);
                hprd_invivo.nodes.put(new_id, new Node(new_id, "protein", "id", "hprd", new_ids));
                n1_invivo = hprd_invivo.nodes.get(new_id);
            } else {
                n1 = hprd.nodes.get(new_id);
                n1_invivo = hprd_invivo.nodes.get(new_id);
            }
            new_id = "p" + conv_map.get(ss[5]);
            if (!hprd.nodes.containsKey(new_id)) {
                new_ids = new ArrayList();
                new_ids.add(new String[]{ss[5]});
                hprd.nodes.put(new_id, new Node(new_id, "protein", "id", "hprd", new_ids));
                n2 = hprd.nodes.get(new_id);
                hprd_invivo.nodes.put(new_id, new Node(new_id, "protein", "id", "hprd", new_ids));
                n2_invivo = hprd.nodes.get(new_id);
            } else {
                n2 = hprd.nodes.get(new_id);
                n2_invivo = hprd_invivo.nodes.get(new_id);
            }
            db_entr = new ArrayList();
            db_entr.add(s);
            hprd.interactions.put("hprd:" + ss[0], new Interaction("hprd:" + ss[0], other_ids, n1, n2, "PPI", "hprd", db_entr, "-", "indirect"));
            if (invivo) {
                hprd_invivo.interactions.put("hprd:" + ss[0], new Interaction("hprd:" + ss[0], other_ids, n1_invivo, n2_invivo, "PPI", "hprd", db_entr, "-", "indirect"));
            }
            invivo = false;
        }
        rd.close();
        System.out.println("hprd.nodes.size " + hprd.nodes.size());
        System.out.println("hprd.interactions.size " + hprd.interactions.size());
        System.out.println("hprd_invivo.nodes.size " + hprd_invivo.nodes.size());
        System.out.println("hprd_invivo.interactions.size " + hprd_invivo.interactions.size());
    }

    /**
     * Load HIPPIE database
     *
     * @param high load all or only high confidence interactions
     */
    public void loadHIPPIE(String dbpath_int_high, String dbpath_int_all, String dbpath_conv, boolean high) throws FileNotFoundException, IOException {
        //entities FoldersPaths.hp_int_high FoldersPaths.hp_int FoldersPaths.hp_conv
        System.out.println("+++++++++++++Hippie++++++++++++");
        String prefix = "";
        List<String> db_entr;
        List<String> other_ids = new ArrayList(); //interactions , сначала пустой
        other_ids.add("-");
        String s;
        String[] ss;
        List<String[]> new_ids;
        String new_id, type;
        BufferedReader rd;
        if (high) {
            rd = new BufferedReader(new FileReader(dbpath_int_high ));
            prefix = "";
        } else {
            rd = new BufferedReader(new FileReader(dbpath_int_all ));
            prefix = "all";
        }
        BufferedReader conv = new BufferedReader(new FileReader(dbpath_conv ));
        Map<String, String> conv_map = new HashMap();
        Node n1, n2, n1_invivo, n2_invivo;
        int i = 0;

        while ((s = conv.readLine()) != null) {
            ss = s.split("\t");
            if (!s.contains("~withdrawn")) {
                conv_map.put(ss[ss.length - 1], ss[0]);
            }
        }
        conv.close();

        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            if (!conv_map.containsKey(ss[1]) || !conv_map.containsKey(ss[3])) {
                continue;
            }
            i++;
            new_id = "p" + conv_map.get(ss[1]);
            if (!hippie.nodes.containsKey(new_id)) {
                new_ids = new ArrayList();
                new_ids.add(new String[]{ss[1]});
                hippie.nodes.put(new_id, new Node(new_id, "protein", "id", "hp", new_ids));
                n1 = hippie.nodes.get(new_id);
            } else {
                n1 = hippie.nodes.get(new_id);
            }
            new_id = "p" + conv_map.get(ss[3]);
            if (!hippie.nodes.containsKey(new_id)) {
                new_ids = new ArrayList();
                new_ids.add(new String[]{ss[3]});
                hippie.nodes.put(new_id, new Node(new_id, "protein", "id", "hp", new_ids));
                n2 = hippie.nodes.get(new_id);
            } else {
                n2 = hippie.nodes.get(new_id);
            }

            db_entr = new ArrayList();
            db_entr.add(s);
            hippie.interactions.put("hp:" + i, new Interaction("hp:" + i, other_ids, n1, n2, "PPI", "hp_" + ss[5], db_entr, ss[4], "indirect"));
            if (Double.valueOf(ss[4]) >= 0.73) {
                //System.out.println(Double.valueOf(ss[4]) +"  "+ ss[4]);
                hippie_high.interactions.put("hp:" + i, new Interaction("hp:" + i, other_ids, n1, n2, "PPI", "hp_" + ss[5], db_entr, ss[4], "indirect"));
                hippie_high.nodes.put(n1.id, n1);
                hippie_high.nodes.put(n2.id, n2);
            }
            if (Double.valueOf(ss[4]) > 0.63) {
                //System.out.println(Double.valueOf(ss[4]) +"  "+ ss[4]);
                hippie_medium_m.interactions.put("hp:" + i, new Interaction("hp:" + i, other_ids, n1, n2, "PPI", "hp_" + ss[5], db_entr, ss[4], "indirect"));
                hippie_medium_m.nodes.put(n1.id, n1);
                hippie_medium_m.nodes.put(n2.id, n2);
            }
            if (Double.valueOf(ss[4]) >= 0.63) {
                //System.out.println(Double.valueOf(ss[4]) +"  "+ ss[4]);
                hippie_medium_meq.interactions.put("hp:" + i, new Interaction("hp:" + i, other_ids, n1, n2, "PPI", "hp_" + ss[5], db_entr, ss[4], "indirect"));
                hippie_medium_meq.nodes.put(n1.id, n1);
                hippie_medium_meq.nodes.put(n2.id, n2);
            }
            if (Double.valueOf(ss[4]) > 0.0) {
                //System.out.println(Double.valueOf(ss[4]) +"  "+ ss[4]);
                hippie_medium_m0.interactions.put("hp:" + i, new Interaction("hp:" + i, other_ids, n1, n2, "PPI", "hp_" + ss[5], db_entr, ss[4], "indirect"));
                hippie_medium_m0.nodes.put(n1.id, n1);
                hippie_medium_m0.nodes.put(n2.id, n2);
            }
            if (s.contains("in vivo")) {
                hippie_invivo.interactions.put("hp:" + i, new Interaction("hp:" + i, other_ids, n1, n2, "PPI", "hp_" + ss[5], db_entr, ss[4], "indirect"));
                hippie_invivo.nodes.put(n1.id, n1);
                hippie_invivo.nodes.put(n2.id, n2);
            }

// hippie.interactions.get("hp:" + i).print();
        }
        rd.close();
        System.out.println("hippie.nodes.size " + hippie.nodes.size());
        System.out.println("hippie.interactions.size " + hippie.interactions.size());
        System.out.println("hippie_high.nodes.size " + hippie_high.nodes.size());
        System.out.println("hippie_high.interactions.size " + hippie_high.interactions.size());
        System.out.println("hippie_medium_m.nodes.size " + hippie_medium_m.nodes.size());
        System.out.println("hippie_medium_m.interactions.size " + hippie_medium_m.interactions.size());
        System.out.println("hippie_medium_meq.nodes.size " + hippie_medium_meq.nodes.size());
        System.out.println("hippie_medium_meq.interactions.size " + hippie_medium_meq.interactions.size());
        System.out.println("hippie_medium_m0.nodes.size " + hippie_medium_m0.nodes.size());
        System.out.println("hippie_medium_m0.interactions.size " + hippie_medium_m0.interactions.size());
        System.out.println("hippie_invivo.size " + hippie_invivo.nodes.size());
        System.out.println("hippie_invivo.size " + hippie_invivo.interactions.size());

    }

    /**
     * Load Signor database
     */
    public void loadSignor(String dbpath_int, String dbpath_conv) throws FileNotFoundException, IOException {
        //entities FoldersPaths.sg_int  FoldersPaths.sg_conv_uniprot_hgnc
        System.out.println("+++++++++++++Signor++++++++++++");
        List<String> db_entr;
        List<String> other_ids = new ArrayList(); //interactions , сначала пустой
        other_ids.add("-");
        String s;
        String[] ss;
        List<String[]> new_ids;
        String new_id, type;
        BufferedReader rd = new BufferedReader(new FileReader(dbpath_int));
        BufferedReader conv_uh = new BufferedReader(new FileReader(dbpath_conv));
        Map<String, String> conv_map_uh = new HashMap();
        Node n1, n2;
        int i = 0;
        int j = 0;
        int k = 0;
        String int_type, prefix, ent_type;
        String id2;
        boolean expr = false;
        while ((s = conv_uh.readLine()) != null) {
            ss = s.split("\t");
            conv_map_uh.put(ss[0], ss[1]);
        }
        conv_uh.close();
        rd.readLine();
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            expr = false;
            i++;
            if (s.contains("STIMULUS")) {
                continue;
            }
            int_type = ss[8] + "_" + ss[9];
            if (int_type.contains("up-regulates quantity by expression") || int_type.contains("transcriptional regulation") || int_type.contains("down-regulates quantity by repression") || int_type.contains("transcriptional repression")) {
                int_type = ss[8] + "_" + ss[9];
                expr = true;
            } else {
                int_type = "PPI_" + ss[8] + "_" + ss[9];
            }

            if (ss[1].equals("PROTEIN")) {
                if (conv_map_uh.containsKey(ss[2])) {
                    new_id = "p" + hugo_by_id.get(conv_map_uh.get(ss[2]))[0];

                } else {
                    new_id = "p" + ss[2];
                }
                if (!signor.nodes.containsKey(new_id)) {
                    new_ids = new ArrayList();
                    new_ids.add(new String[]{conv_map_uh.get(ss[2])});
                    signor.nodes.put(new_id, new Node(new_id, "protein", "id", "sg", new_ids));
                    n1 = signor.nodes.get(new_id);
                } else {
                    n1 = signor.nodes.get(new_id);
                }
            } else {
                new_id = ss[2];
                if (!signor.nodes.containsKey(new_id)) {
                    new_ids = new ArrayList();
                    new_ids.add(new String[]{ss[2]});
                    signor.nodes.put(new_id, new Node(new_id, ss[1], "id", "sg", new_ids));
                    n1 = signor.nodes.get(new_id);
                } else {
                    n1 = signor.nodes.get(new_id);
                }
            }

            if (ss[5].equals("PROTEIN")) {
                if (expr) {
                    prefix = "";
                    ent_type = "gene";
                } else {
                    prefix = "p";
                    ent_type = "protein";
                }
                if (conv_map_uh.containsKey(ss[6])) {
                    new_id = prefix + hugo_by_id.get(conv_map_uh.get(ss[6]))[0];
                } else {
                    new_id = prefix + ss[6];
                }
                if (!signor.nodes.containsKey(new_id)) {
                    new_ids = new ArrayList();
                    new_ids.add(new String[]{conv_map_uh.get(ss[6])});
                    signor.nodes.put(new_id, new Node(new_id, ent_type, "id", "sg", new_ids));
                    n2 = signor.nodes.get(new_id);
                } else {
                    n2 = signor.nodes.get(new_id);
                }
            } else {
                new_id = ss[6];
                if (expr) {
                    System.out.println(s);
                    //break;
                }
                if (!signor.nodes.containsKey(new_id)) {
                    new_ids = new ArrayList();
                    new_ids.add(new String[]{ss[6]});
                    signor.nodes.put(new_id, new Node(new_id, ss[5], "id", "sg", new_ids));
                    n2 = signor.nodes.get(new_id);
                } else {
                    n2 = signor.nodes.get(new_id);
                }
            }

            //System.out.print(hugo_by_id.get(n1.id.substring(1))[1]);
            //System.out.print("\t" + hugo_by_id.get(n2.id.substring(1))[1] + "\t");
            // if (ss[1].equals("PROTEINFAMILY") || ss[1].equals("SMALLMOLECULE") || ss[1].equals("CHEMICAL") || ss[1].equals("COMPLEX") || ss[1].equals("PHENOTYPE")) {
            // if (ss[5].equals("PROTEINFAMILY") || ss[5].equals("SMALLMOLECULE") || ss[5].equals("CHEMICAL") || ss[5].equals("COMPLEX") || ss[5].equals("PHENOTYPE")) {
            db_entr = new ArrayList();
            db_entr.add(s);
            signor.interactions.put("sg:" + i, new Interaction("sg:" + i, other_ids, n1, n2, int_type, "sg", db_entr, "1", "direct"));
            //signor.interactions.get("sg:" + i).print();
        }
        rd.close();
        System.out.println("signor.nodes.size " + signor.nodes.size());
        System.out.println("signor.interactions.size " + signor.interactions.size());
    }

    /**
     * Load Signalink database
     */
    public void loadSignalink(String dbpath_int, String dbpath_conv) throws FileNotFoundException, IOException {
        //entities FoldersPaths.sl_int FoldersPaths.sl_conv_uniprot_hgnc
        System.out.println("+++++++++++++Signalink++++++++++++");
        List<String> db_entr;
        List<String> other_ids = new ArrayList(); //interactions , сначала пустой
        other_ids.add("-");
        String s;
        String[] ss;
        List<String[]> new_ids;
        String new_id, direct;
        BufferedReader rd = new BufferedReader(new FileReader(dbpath_int));
        BufferedReader conv_uh = new BufferedReader(new FileReader(dbpath_conv));
        Map<String, String> conv_map_uh = new HashMap();
        Node n1, n2;
        int i = 0;
        int j = 0;
        int k = 0;
        String id1;
        String id2;
        while ((s = conv_uh.readLine()) != null) {
            ss = s.split("\t");
            conv_map_uh.put(ss[0], ss[1]);
        }
        conv_uh.close();
        rd.readLine();
        String id_type;
        while ((s = rd.readLine()) != null) {
            id_type = "";
            ss = s.split(";");
            i++;
            if (!ss[3].contains("H. sapiens") || !ss[9].contains("H. sapiens") || (s.contains("ELM based prediction"))) {
                continue;
            }

            if (ss[13].contains("PPI directed") || ss[13].contains("PPI predicted as directed") || ss[13].contains("Transcriptional directed") || ss[13].contains("Post-transcriptional directed")) {
                direct = "direct";
            } else {
                direct = "indirect";
            }

            if (ss[13].contains("PPI directed") || ss[13].contains("PPI predicted as directed") || ss[13].contains("PPI undirected")) {

                if (conv_map_uh.containsKey(ss[1])) {
                    new_id = "p" + hugo_by_id.get(conv_map_uh.get(ss[1]))[0];
                } else {
                    new_id = "p" + ss[1];
                }
                if (!signalink.nodes.containsKey(new_id)) {
                    new_ids = new ArrayList();
                    new_ids.add(new String[]{conv_map_uh.get(ss[1])});
                    signalink.nodes.put(new_id, new Node(new_id, "protein", "id", "signalink", new_ids));
                }
                n1 = signalink.nodes.get(new_id);
                if (conv_map_uh.containsKey(ss[7])) {
                    new_id = "p" + hugo_by_id.get(conv_map_uh.get(ss[7]))[0];
                } else {
                    new_id = "p" + ss[7];
                }
                if (!signalink.nodes.containsKey(new_id)) {
                    new_ids = new ArrayList();
                    new_ids.add(new String[]{conv_map_uh.get(ss[7])});
                    signalink.nodes.put(new_id, new Node(new_id, "protein", "id", "signalink", new_ids));
                }
                n2 = signalink.nodes.get(new_id);
            } else if (ss[13].contains("Transcriptional directed")) {
                if (conv_map_uh.containsKey(ss[1])) {
                    new_id = "p" + hugo_by_id.get(conv_map_uh.get(ss[1]))[0];
                } else {
                    new_id = "p" + ss[1];
                }
                if (!signalink.nodes.containsKey(new_id)) {
                    new_ids = new ArrayList();
                    new_ids.add(new String[]{conv_map_uh.get(ss[1])});
                    signalink.nodes.put(new_id, new Node(new_id, "protein", "id", "signalink", new_ids));
                }
                n1 = signalink.nodes.get(new_id);
                if (conv_map_uh.containsKey(ss[7])) {
                    new_id = hugo_by_id.get(conv_map_uh.get(ss[7]))[0];
                } else {
                    new_id = ss[7];
                }
                if (!signalink.nodes.containsKey(new_id)) {
                    new_ids = new ArrayList();
                    new_ids.add(new String[]{conv_map_uh.get(ss[7])});
                    signalink.nodes.put(new_id, new Node(new_id, "gene", "id", "signalink", new_ids));

                }
                n2 = signalink.nodes.get(new_id);
                id_type = "_tf";
            } else {
                new_id = ss[0];
                if (!signalink.nodes.containsKey(new_id)) {
                    new_ids = new ArrayList();
                    new_ids.add(new String[]{"-"});
                    signalink.nodes.put(new_id, new Node(new_id, "miRNA", "id", "signalink", new_ids));
                }
                n1 = signalink.nodes.get(new_id);
                id_type = "_mr";
                if (conv_map_uh.containsKey(ss[7])) {
                    new_id = hugo_by_id.get(conv_map_uh.get(ss[7]))[0];
                } else {
                    new_id = ss[7];
                }
                if (!signalink.nodes.containsKey(new_id)) {
                    new_ids = new ArrayList();
                    new_ids.add(new String[]{conv_map_uh.get(ss[7])});
                    signalink.nodes.put(new_id, new Node(new_id, "gene", "id", "signalink", new_ids));

                }
                n2 = signalink.nodes.get(new_id);

            }

            //System.out.print(hugo_by_id.get(n1.id.substring(1))[1]);
            //System.out.print("\t" + hugo_by_id.get(n2.id.substring(1))[1] + "\t");
            // if (ss[1].equals("PROTEINFAMILY") || ss[1].equals("SMALLMOLECULE") || ss[1].equals("CHEMICAL") || ss[1].equals("COMPLEX") || ss[1].equals("PHENOTYPE")) {
            // if (ss[5].equals("PROTEINFAMILY") || ss[5].equals("SMALLMOLECULE") || ss[5].equals("CHEMICAL") || ss[5].equals("COMPLEX") || ss[5].equals("PHENOTYPE")) {
            db_entr = new ArrayList();
            db_entr.add(s);
            signalink.interactions.put("sl:" + i + id_type, new Interaction("sl:" + i, other_ids, n1, n2, ss[13], "sl_" + ss[17], db_entr, "1", direct));
            // signalink.interactions.get("sl:" + i).print();
        }
        rd.close();
        System.out.println("signalink.nodes.size " + signalink.nodes.size());
        System.out.println("signalink.interactions.size " + signalink.interactions.size());
    }

    /**
     * Load tFacts database
     */
    public void loadtFacts(String dbpath_int) throws FileNotFoundException, IOException {
        // FoldersPaths.tfactsf
        System.out.println("+++++++++++++tFacts++++++++++++");
        String s;
        List<String> db_entr;
        List<String> other_ids = new ArrayList(); //interactions , сначала пустой
        other_ids.add("-");
        String[] ss;// = new String[5];
        List<String[]> new_ids;
        String new_id;// = new String();
        Node n1, n2;
        Interaction int1;
        String id_type;
        String db_flag;
        int i = 0;
        int j = 0;
        int k = 0;
        //entities  & interactions       
        BufferedReader rd = new BufferedReader(new FileReader(dbpath_int));
        s = rd.readLine();
        while ((s = rd.readLine()) != null) {
            i++;
            ss = s.split("\t");
            if (ss[3].toLowerCase().contains("human") || ss[3].toLowerCase().contains("homo sapiens")) {
                j++;
                //Target gene ss[0]                 
                if (hugo.containsKey(ss[0])) {
                    new_id = hugo.get(ss[0])[0];
                    id_type = "hugo";
                    db_flag = "1000007";
                } else {
                    if (tfacts_to_hugo.containsKey(ss[0])) {
                        new_id = tfacts_to_hugo.get(ss[0]);
                        id_type = "hugo";
                        db_flag = "1000007";
                    } else {
                        new_id = ss[0];
                        id_type = "tfacts";
                        db_flag = "0000007";
                    }
                }
                if (!tfacts.nodes.containsKey(new_id)) {
                    new_ids = new ArrayList();
                    new_ids.add(new String[]{ss[0]});
                    //new_ids.add(hugo.get(ss[0]));                    
                    new_ids.add(ss);
                    n1 = new Node(new_id, "gene", id_type, "tf", new_ids);
                    tfacts.nodes.put(new_id, n1);
                } else {
                    n1 = tfacts.nodes.get(new_id);
                }
                //transcription factor ss[1]
                if (hugo.containsKey(ss[1])) {
                    new_id = "p" + hugo.get(ss[1])[0];
                    id_type = "hugo";
                    db_flag = "1000007";
                    //System.out.println("=1= " + new_id);
                } else {
                    if (tfacts_to_hugo.containsKey(ss[1])) {
                        new_id = "p" + tfacts_to_hugo.get(ss[1]);
                        id_type = "hugo";
                        db_flag = "1000007";
                        // System.out.println("=2= " + new_id);
                    } else {
                        new_id = ss[1];
                        id_type = "tfacts";
                        db_flag = "0000007";
                        // System.out.println("=3= " + new_id);
                    }
                }
                if (!tfacts.nodes.containsKey(new_id)) {
                    // System.out.println("=3= " + new_id);
                    new_ids = new ArrayList();
                    // new_ids.add(hugo.get(ss[1]));
                    new_ids.add(new String[]{ss[1]});
                    new_ids.add(ss);
                    n2 = new Node(new_id, "protein", id_type, "tf", new_ids);
                    tfacts.nodes.put(new_id, n2);
                } else {
                    n2 = tfacts.nodes.get(new_id);
                    //System.out.println("=2= " + new_id + "\t" + n2.id);
                }
                db_entr = new ArrayList();
                db_entr.add(s);
                int1 = new Interaction("tfacts:" + i, other_ids, n2, n1, "TF-gene ", "tfacts", db_entr, ss[4], "direct");
                tfacts.interactions.put("tfacts:" + i, int1);
                //System.out.println("tfacts:" + i  + "\t" + n2.id + "\t" + n1.id);
            }
        }
        rd.close();
        System.out.println("i=all " + i + " j=human " + j + " k= " + k);
        System.out.println("tfacts_int " + tfacts.interactions.size() + " tfacts_ent " + tfacts.nodes.size());
    }

    /**
     * Load KEGG database
     */
    public void loadKEGG(String dbpath_int, String dbpath_ent) throws FileNotFoundException, IOException {
        // FoldersPaths.kegg_ent_comp FoldersPaths.kegg_int
        System.out.println("+++++++++++++KEGG++++++++++++");
        List<String> db_entr;
        List<String> other_ids = new ArrayList(); //interactions , сначала пустой
        Network kegg_nodes = new Network();
        other_ids.add("-");
        String s;
        String[] ss = new String[2];
        List<String[]> new_ids;
        String new_id = new String();
        Node n1, n2;
        int i = 0;
        int j = 0;
        int k = 0;
        //entities compounds and glycans only
        BufferedReader rd = new BufferedReader(new FileReader(dbpath_ent));
        while ((s = rd.readLine()) != null) {
            i++;
            new_ids = new ArrayList();
            ss = s.split("\t");
            new_id = ss[0];
            new_ids.add(ss);
            //new_ids.add(ss);
            kegg_nodes.nodes.put(new_id, new Node(new_id, "compound", "id", "kegg", new_ids));
        }
        rd.close();
        rd = new BufferedReader(new FileReader(FoldersPaths.kegg_ent_glyc));
        while ((s = rd.readLine()) != null) {
            j++;
            new_ids = new ArrayList();
            ss = s.split("\t");
            new_id = ss[0];
            //new_ids.add(null);
            new_ids.add(ss);
            kegg_nodes.nodes.put(new_id, new Node(new_id, "glycan", "id", "kegg", new_ids));
        }
        rd.close();
        //interactions and genes
        ss = new String[11];
        String new_id1, new_id2;
        rd = new BufferedReader(new FileReader(dbpath_int));
        while ((s = rd.readLine()) != null) {
            new_ids = new ArrayList();
            ss = s.split("\t");
            new_id1 = ss[4];
            new_id2 = ss[5];
            if (ss[4].contains("HGNC:")) {
                new_id1 = "p" + ss[4];
            }
            if (ss[4].contains("hsa:")) {
                new_id1 = ss[4];
                if (!kegg_nodes.nodes.containsKey(ss[4])) {
                    //new_ids.add(null);
                    new_ids.add(ss);
                    kegg_nodes.nodes.put(ss[4], new Node(ss[4], "protein", "KEGG_entrez", "kegg", new_ids));
                }
            }
            if (ss[4].contains("HGNC:") && (!kegg_nodes.nodes.containsKey("p" + ss[4]))) {
                k++;
                new_ids.add(hugo_by_id.get(ss[4]));
                kegg_nodes.nodes.put("p" + ss[4], new Node("p" + ss[4], "protein", "hugo", "kegg", new_ids));
            }
            new_ids = new ArrayList();

            if (ss[5].contains("HGNC:")) {
                new_id2 = "p" + ss[5];
            }
            if (ss[5].contains("hsa:")) {
                new_id2 = ss[5];
                if (!kegg_nodes.nodes.containsKey(ss[5])) {
                    //new_ids.add(null);
                    new_ids.add(ss);
                    kegg_nodes.nodes.put(ss[5], new Node(ss[5], "protein", "KEGG_entrez", "kegg", new_ids));
                }
            }
            if (ss[5].contains("HGNC:") && (!kegg_nodes.nodes.containsKey("p" + ss[5]))) {
                k++;
                new_ids.add(hugo_by_id.get(ss[5]));
                kegg_nodes.nodes.put("p" + ss[5], new Node("p" + ss[5], "protein", "hugo", "kegg", new_ids));
            }
            if (!ss[4].contains("HGNC:") && !ss[4].contains("cpd:") && !ss[4].contains("gl:") && !ss[5].contains("HGNC:") && !ss[5].contains("cpd:") && !ss[5].contains("gl:")) {
                System.out.println("problem with ids in kegg " + "s");
            }
            db_entr = new ArrayList();
            db_entr.add(s);
            try {
                String q = kegg_nodes.nodes.get(new_id2).id;
            } catch (NullPointerException e) {
                System.out.println("zzzz " + new_id2 + "  " + s);
                // break;
            }
            kegg.interactions.put(ss[0], new Interaction(ss[0], other_ids, kegg_nodes.nodes.get(new_id1), kegg_nodes.nodes.get(new_id2), ss[3], "kegg", db_entr, "kegg", ss[10]));
            kegg.nodes.put(new_id1, kegg_nodes.nodes.get(new_id1));
            kegg.nodes.put(new_id2, kegg_nodes.nodes.get(new_id2));
        }
        rd.close();
        // kegg.nodes=kegg_nodes.nodes;
        System.out.println("i=compounds " + i + " j=glycans " + j + " k=genes " + k);
        System.out.println("kegg_int " + kegg.interactions.size() + " kegg_ent " + kegg.nodes.size());
    }

    /**
     * Load TransMir database
     */
    public void loadTransmir(String dbpath_int, String dbpath_conv) throws FileNotFoundException, IOException {
        //  FoldersPaths.trsmir_crstbl  FoldersPaths.trsmir
        System.out.println("+++++++++++++Transmir++++++++++++");
        String s;
        List<String> db_entr;
        List<String> other_ids = new ArrayList(); //interactions , сначала пустой
        other_ids.add("-");
        String[] ss;// = new String[10];
        List<String[]> new_ids;
        String new_id;// = new String();
        Node n1, n2;
        Interaction int1;
        int ii = 0;
        int j = 0;
        int k = 0;
        String mirna_gene;
        int temp;
        //hugo_mirna
        List<String> mirna_gene_names;
        Map<String, List<String>> mirna_hugo_crtbl = new HashMap();
        List<String> tmp_list;
        BufferedReader rd = new BufferedReader(new FileReader(dbpath_conv));
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            if (mirna_hugo_crtbl.containsKey(ss[0])) {
                tmp_list = mirna_hugo_crtbl.get(ss[0]);
                tmp_list.add(ss[1]);
                mirna_hugo_crtbl.put(ss[0], tmp_list);
            } else {
                tmp_list = new ArrayList();
                tmp_list.add(ss[1]);
                mirna_hugo_crtbl.put(ss[0], tmp_list);
            }
        }
        rd.close();
        //entities  & interactions       
        rd = new BufferedReader(new FileReader(dbpath_int));
        s = rd.readLine();
        while ((s = rd.readLine()) != null) {
            ii++;
            ss = s.split("\t");
            if (ss[9].toLowerCase().equals("human")) {
                j++;
                //TF  gene
                if (!ss[1].equals("na")) {
                    try {
                        new_id = "p" + hugo_entrez.get(ss[1])[0];
                    } catch (NullPointerException e) {
                        System.out.println(" ------- No HUGO name for tf: " + ss[0]);
                        continue;
                    }
                    new_ids = new ArrayList();
                    new_ids.add(ss);
                    new_ids.add(hugo_entrez.get(ss[1]));

                    n1 = new Node(new_id, "protein", "hugo", "tm", new_ids); //тип был TF                    
                    if (transmir.nodes.containsKey(new_id)) {
                        transmir.nodes.get(new_id).ids.add(ss);
                    } else {
                        transmir.nodes.put(new_id, n1);
                    }
                    //mirna                    
                    // ----------------11/2015 new
                    if (ss[3].contains("hsa-")) {
                        mirna_gene = ss[3].substring(4);
                    } else {
                        mirna_gene = ss[3];
                    }
                    if (mirna_gene.contains("let")) {
                        mirna_gene = "MIR" + mirna_gene.replaceFirst("-", "").toUpperCase();
                    } else {
                        mirna_gene = mirna_gene.replaceFirst("-", "").toUpperCase();
                    }
                    if (mirna_gene.contains("-")) {
                        temp = mirna_gene.indexOf("-");
                        if (Character.isLetter(mirna_gene.charAt(temp - 1)) && Character.isDigit(mirna_gene.charAt(temp + 1))) {
                            //System.out.println(mirna_gene);
                            mirna_gene = mirna_gene.replaceFirst("-", "");
                        }
                    }
                    mirna_gene_names = new ArrayList();
                    if (hugo.containsKey(mirna_gene)) {
                        mirna_gene_names.add(hugo.get(mirna_gene)[0]);
                        // System.out.println(ss[3] + "\t" + hugo.get(mirna_gene)[0]);
                    } else {
                        if (mirna_hugo_crtbl.containsKey(ss[3])) {
                            mirna_gene_names = mirna_hugo_crtbl.get(ss[3]);
                        } else {
                            System.out.println("------- No HUGO gene for the miRNA: " + ss[3]);
                        }
                    }
                    // -------------new
                  /*  for (String t : mirna_gene_names) {
                     i++;
                     new_ids = new ArrayList();
                     new_ids.add(null);
                     new_ids.add(ss);
                     n2 = new Node(t, "gene", "transmir", "10005", new_ids);
                     if (transmir.nodes.containsKey("hsa-" + ss[3])) {
                     transmir.nodes.get(t).ids.add(ss);
                     } else {
                     transmir.nodes.put(t, n2);
                     }
                     //interaction
                     db_entr = new ArrayList();
                     db_entr.add(ss);
                     int1 = new Interaction("transmir:" + i, other_ids, n1, n2, "TF-miRNA gene" + ss[7], "transmir", db_entr, "", ss[8]);
                     transmir.interactions.put("transmir:" + i, int1);
                     //System.out.println ("transmir:" + i + "\t"+ss[0]+"\t"+ n1.id + "\t" +hugo_by_id.get(n1.id.substring(1))[1] + "\t"+ss[3]+ "\t" + n2.id+ "\t" +hugo_by_id.get(n2.id)[1]);                        
                     } */
                    //i++;
                    new_ids = new ArrayList();
                    // new_ids.add(null);
                    new_ids.add(ss);
                    n2 = new Node("hsa-" + ss[3].toLowerCase(), "mirna", "transmir", "tm", new_ids);
                    if (transmir.nodes.containsKey("hsa-" + ss[3].toLowerCase())) {
                        transmir.nodes.get("hsa-" + ss[3].toLowerCase()).ids.add(ss);
                    } else {
                        transmir.nodes.put("hsa-" + ss[3].toLowerCase(), n2);
                    }
                    //interaction
                    db_entr = new ArrayList();
                    db_entr.add(s);
                    int1 = new Interaction("transmir:" + ii, other_ids, n1, n2, "TF-miRNA gene" + ss[7], "transmir", db_entr, "", ss[8]);
                    transmir.interactions.put("transmir:" + ii, int1);
                    //System.out.println("transmir:" + ii + "\t" + n1.id + "\t" +n2.id);
                } else {
                    System.out.println(s);
                }
            }
        }
        rd.close();
        System.out.println("i=all " + ii + " j=human " + j + " k= " + k);
        System.out.println("transmir_int " + transmir.interactions.size() + " transmir_ent " + transmir.nodes.size());
    }

    /**
     * Load mitTarBase database
     */
    public void loadmirTarBase(String dbpath_int) throws FileNotFoundException, IOException {
        // FoldersPaths.mirtarb
        System.out.println("+++++++++++++mirTarBase++++++++++++");
        String s;
        String[] ss = new String[9];
        List<String[]> new_ids;
        List<String> db_entr;
        List<String> other_ids = new ArrayList(); //interactions , сначала пустой
        other_ids.add("-");
        String new_id = new String();
        Node n1, n2;
        Interaction int1;
        int i = 0;
        int j = 0;
        int k = 0;

        //entities  & interactions       
        BufferedReader rd = new BufferedReader(new FileReader(dbpath_int));
        s = rd.readLine();
        while ((s = rd.readLine()) != null) {
            i++;
            ss = s.split("\t");
            if (ss[5].toLowerCase().equals("homo sapiens")) {
                j++;
                //mirna
                if (ss[1].endsWith("-3p") || ss[1].endsWith("-5p")) {
                    new_id = ss[1].substring(0, ss[1].length() - 3).toLowerCase();
                    // System.out.println(new_id);
                } else if (ss[1].endsWith("*")) {
                    new_id = ss[1].substring(0, ss[1].length() - 1).toLowerCase();
                    // System.out.println(new_id);
                } else {
                    new_id = ss[1].toLowerCase();
                }
                new_ids = new ArrayList();
                // new_ids.add(null);
                new_ids.add(ss);
                n1 = new Node(new_id, "miRNA", "mirtarbase", "mtb", new_ids);
                if (mirtarbase.nodes.containsKey(new_id)) {
                    mirtarbase.nodes.get(new_id).ids.add(ss);
                } else {
                    mirtarbase.nodes.put(new_id, n1);
                }
                String id_type;
                //gene    
                if (ss[4] != null) {
                    new_ids = new ArrayList();
                    try {
                        new_id = hugo_entrez.get(ss[4])[0];
                        id_type = "hugo";
                        new_ids.add(hugo_entrez.get(ss[4]));
                    } catch (NullPointerException e) {
                        new_id = "entrz:" + ss[4];
                        id_type = "entrez";
                        // new_ids.add(null);
                        new_ids.add(ss);
                    }
                    new_ids.add(ss);
                    n2 = new Node(new_id, "gene", id_type, "mtb", new_ids);
                    if (mirtarbase.nodes.containsKey(new_id)) {
                        mirtarbase.nodes.get(new_id).ids.add(ss);
                    } else {
                        mirtarbase.nodes.put(new_id, n2);
                    }
                    db_entr = new ArrayList();
                    db_entr.add(s);
                    int1 = new Interaction(ss[0], other_ids, n1, n2, "miRNA-gene ", "mirtarbase", db_entr, ss[6], ss[8]);
                    if (mirtarbase.interactions.containsKey(ss[0])) {
                    }
                    mirtarbase.interactions.put(ss[0], int1);
                } else {
                    System.out.println(s);
                }
            }
        }
        rd.close();
        System.out.println("i=all " + i + " j=human " + j + " k= " + k);
        System.out.println("mirtarbase_int " + mirtarbase.interactions.size() + " mirtarbase_ent " + mirtarbase.nodes.size());
    }
}


