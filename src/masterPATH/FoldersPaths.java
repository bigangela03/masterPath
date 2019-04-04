package masterPATH;

/**
 * FolderPaths class contains system paths to databases files, files with pathways, paths, nodes, etc
 * 
 * @author Natalia Rubanova
 */
public class FoldersPaths {

   
    static final String dbpath = "";
    static final String syspath = "";
    static final String sep = "";
     


    static final String hprd_folder = "" + sep;
    static final String transpath_folder = "" + sep;
    static final String kegg_folder = "" + sep;
    static final String transmir_folder = "" + sep;
    static final String mirtarbae_folder = "" + sep;
    static final String tfacts_folder = "" + sep;
    static final String hugo_folder = "" + sep;
    static String[] hugo_path = {dbpath + hugo_folder + "hugo_all_17Jun.txt"};

    static final String hp_folder = "" + sep;
    static final String sg_folder = "" + sep;
    static final String sl_folder = "" + sep;

    static String hugo_with_entrez = dbpath + hugo_folder + "HUGO_with_entrezid.txt";

    static final String hp_int = dbpath + hp_folder + "hippie_current.txt";
    static final String hp_int_high = dbpath + hp_folder + "hippie_high.txt";
    static final String hp_conv = dbpath + hp_folder + "hippie_hgnc_to_entrz.txt";

    static final String sg_int = dbpath + sg_folder + "all_data_edt.tsv";
    static final String sg_conv_uniprot_hgnc = dbpath + sg_folder + "signor_uniprot_hgnc.txt";

    static final String sl_int = dbpath + sl_folder + "06162016-signalink-VA5Eui_edt.csv";
    static final String sl_conv_uniprot_hgnc = dbpath + sl_folder + "signalink_uniprot_hgnc.txt";

    static String hprd_ent = dbpath + hprd_folder + "PROTEIN_NOMENCLATURE_EDT.txt";
    static String hprd_int = dbpath + hprd_folder + "BINARY_PROTEIN_PROTEIN_INTERACTIONS_EDT.txt";
    static String hprd_bioDB = dbpath + hprd_folder + "missed_hprd_to_hugo_bioDB.txt";
    static String hprd_conv = dbpath + hprd_folder + "hprd_id_conv.txt";


    static String kegg_ent_comp = dbpath + kegg_folder + "list_of_all_compounds.txt";
    static String kegg_ent_glyc = dbpath + kegg_folder + "list_of_all_glycans.txt";
    static String kegg_int = dbpath + kegg_folder + "mergedKEGGinteractions.txt";
    static String trsmir = dbpath + transmir_folder + "transmir_v1.2.txt";
    static String trsmir_crstbl = dbpath + transmir_folder + "miRNA_hugo.txt";
    static String mirtarb = dbpath + mirtarbae_folder + "hsa_MTI.txt";
    static String entrz_wo_hugo = dbpath + mirtarbae_folder + "entrez_wo_hugo.txt";
    static String tfactsf = dbpath + tfacts_folder + "tfacts.txt";
    static String tfacts_crosstable = dbpath + tfacts_folder + "tfacts_crosstable.txt";

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////SYSTEMS/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
    static final String ogg1 = "OGG1" + sep;
    static final String inputOG = "input" + sep;
    static final String outputOG = "output" + sep;
    static final String outputFOG = syspath + ogg1 + outputOG;
    static String hlOG = syspath + ogg1 + inputOG + "HitGenes.txt";
    static String fpOG = syspath + ogg1 + inputOG + "FinalPlayers.txt";
    static String fpOG_ppi = syspath + ogg1 + inputOG + "FinalPlayers3.txt";
    static String foundDFOG = syspath + ogg1 + outputOG + "foundDF_28_11.txt";
    static String foundDFOG_overrepr = syspath + ogg1 + outputOG + "foundDF_overrepr.txt";
    static String foundDFOG_mirnas = syspath + ogg1 + outputOG + "foundDFmirnas.txt";
    static String OGrandom_path = syspath + ogg1 + outputOG + "random" + sep;
    static final String ppiOG = "ppi" + sep;
    static String foundDFOG_ppi = syspath + ogg1 + outputOG + ppiOG + "foundDF_18_01.txt";
    static String OGrandom_path_ppi = syspath + ogg1 + outputOG + ppiOG + "random" + sep;
    static String OGrandom_path2_ppi = syspath + ogg1 + outputOG + ppiOG + "random2" + sep;

    static final String mirnas63 = "63mirnas" + sep;
    static final String input63 = "input" + sep;
    static final String output63 = "output" + sep;
    static final String outputF63 = syspath + mirnas63 + output63;
    static String hl63 = syspath + mirnas63 + input63 + "63mirnas_clean.txt";
    static String fp63 = syspath + mirnas63 + input63 + "FinalPlayers.txt";
    static String foundDF63 = syspath + mirnas63 + output63 + "foundDF_18_01.txt";
    static String foundDF63_overrepr = syspath + mirnas63 + output63 + "foundDF_overrepr.txt";
    static String foundDF63_mirnas = syspath + mirnas63 + output63 + "foundDFmirnas.txt";
    static String M63random_path = syspath + mirnas63 + output63 + "random" + sep;

    static final String case2 = "Case2_Control" + sep;
    static final String inputC2 = "input" + sep;
    static final String outputC2 = "output" + sep;
    static final String outputFC2 = syspath + case2 + outputC2;
    static String hlC2 = syspath + case2 + inputC2 + "HitGenes.txt";
    static String fpC2 = syspath + case2 + inputC2 + "FinalPlayers.txt";
    static String exprC2 = syspath + case2 + inputC2 + "expression.txt";
    static String foundDFC2 = syspath + case2 + outputC2 + "foundDF_18_01.txt";
    static String foundDFC2_overrepr = syspath + case2 + outputC2 + "foundDF_overrepr.txt";
    static String foundDFC2_mirnas = syspath + case2 + outputC2 + "foundDFmirnas.txt";
    static String C2random_path = syspath + case2 + outputC2 + "random" + sep;
    static String C2random_path2 = syspath + case2 + outputC2 + "random2" + sep;

    

    static final String stat = "Statistics" + sep;
    static final String inputstatf = "input" + sep;
    static final String outputstatf = "output" + sep;
    static final String outputstat = syspath + stat + outputstatf;

}


    /*
    static final String all_nodes = dbpath + "all_nodes.txt";
    static final String all_interactions = dbpath + "all_interactions.txt";
    static final String pc_folder = "PC" + sep;
    static final String pc_int = dbpath + pc_folder + "PathwayCommons.8.All.BINARY_SIF.hgnc.txt.sif"; */


    /*   static final String dbpath = "//home//nrubanova//Documents//PhD//DB//";
     static final String syspath = "//home//nrubanova//Documents//PhD//Systems//";
     static final String sep = "//"; */

/*
    static String tp_ent = dbpath + transpath_folder + "reaction_semantic_full_entities_EDT.txt";
    static String tp_int = dbpath + transpath_folder + "reaction_semantic_full_connections_EDT.txt";
*/


    // static String[] hugo_path = {dbpath + hugo_folder + "protein-coding_gene.txt", dbpath + hugo_folder + "other.txt", dbpath + hugo_folder + "non-coding_RNA.txt", dbpath + hugo_folder + "pseudogene.txt", dbpath + hugo_folder + "phenotype.txt"};

/*

    static final String hprd_folder = "HPRD" + sep;
    static final String transpath_folder = "Transpath" + sep;
    static final String kegg_folder = "KEGG_03_15" + sep;
    static final String transmir_folder = "transmir" + sep;
    static final String mirtarbae_folder = "mirtarbase" + sep;
    static final String tfacts_folder = "tfacts" + sep;
    static final String hugo_folder = "GeneNames_HUGO" + sep;
    static String[] hugo_path = {dbpath + hugo_folder + "hugo_all_17Jun.txt"};

    static final String hp_folder = "Hippie" + sep;
    static final String sg_folder = "Signor" + sep;
    static final String sl_folder = "SignaLink" + sep;

    static String hugo_with_entrez = dbpath + hugo_folder + "HUGO_with_entrezid.txt";

    static final String hp_int = dbpath + hp_folder + "hippie_current.txt";
    static final String hp_int_high = dbpath + hp_folder + "hippie_high.txt";
    static final String hp_conv = dbpath + hp_folder + "hippie_hgnc_to_entrz.txt";

    static final String sg_int = dbpath + sg_folder + "all_data_edt.tsv";
    static final String sg_conv_uniprot_hgnc = dbpath + sg_folder + "signor_uniprot_hgnc.txt";

    static final String sl_int = dbpath + sl_folder + "06162016-signalink-VA5Eui_edt.csv";
    static final String sl_conv_uniprot_hgnc = dbpath + sl_folder + "signalink_uniprot_hgnc.txt";

    static String hprd_ent = dbpath + hprd_folder + "PROTEIN_NOMENCLATURE_EDT.txt";
    static String hprd_int = dbpath + hprd_folder + "BINARY_PROTEIN_PROTEIN_INTERACTIONS_EDT.txt";
    static String hprd_bioDB = dbpath + hprd_folder + "missed_hprd_to_hugo_bioDB.txt";
    static String hprd_conv = dbpath + hprd_folder + "hprd_id_conv.txt";


    static String kegg_ent_comp = dbpath + kegg_folder + "list_of_all_compounds.txt";
    static String kegg_ent_glyc = dbpath + kegg_folder + "list_of_all_glycans.txt";
    static String kegg_int = dbpath + kegg_folder + "mergedKEGGinteractions.txt";
    static String trsmir = dbpath + transmir_folder + "transmir_v1.2.txt";
    static String trsmir_crstbl = dbpath + transmir_folder + "miRNA_hugo.txt";
    static String mirtarb = dbpath + mirtarbae_folder + "hsa_MTI.txt";
    static String entrz_wo_hugo = dbpath + mirtarbae_folder + "entrez_wo_hugo.txt";
    static String tfactsf = dbpath + tfacts_folder + "tfacts.txt";
    static String tfacts_crosstable = dbpath + tfacts_folder + "tfacts_crosstable.txt";


    //static final String dbpath = "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\DB\\";
    //static final String syspath = "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\";
    //static final String sep = "\\";

*/