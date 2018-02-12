package masterPATH;

/**
 * FolderPaths class contains system paths to databases files, files with pathways, paths, nodes, etc
 * 
 * @author Natalia Rubanova
 */
public class FoldersPaths {

    static final String dbpath = "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\DB\\";
    static final String syspath = "F:\\Dropbox\\_Работа\\Programms\\DATA\\PhD\\Systems\\";
    static final String sep = "\\";
    /*   static final String dbpath = "//home//nrubanova//Documents//PhD//DB//";
     static final String syspath = "//home//nrubanova//Documents//PhD//Systems//";
     static final String sep = "//"; */

    static final String hprd_folder = "HPRD" + sep;
    static final String transpath_folder = "Transpath" + sep;
    static final String kegg_folder = "KEGG_03_15" + sep;
    static final String transmir_folder = "transmir" + sep;
    static final String mirtarbae_folder = "mirtarbase" + sep;
    static final String tfacts_folder = "tfacts" + sep;
    static final String hugo_folder = "GeneNames_HUGO" + sep;
    // static String[] hugo_path = {dbpath + hugo_folder + "protein-coding_gene.txt", dbpath + hugo_folder + "other.txt", dbpath + hugo_folder + "non-coding_RNA.txt", dbpath + hugo_folder + "pseudogene.txt", dbpath + hugo_folder + "phenotype.txt"};
    static String[] hugo_path = {dbpath + hugo_folder + "hugo_all_17Jun.txt"};
    static final String all_nodes = dbpath + "all_nodes.txt";
    static final String all_interactions = dbpath + "all_interactions.txt";
    static final String pc_folder = "PC" + sep;
    static final String pc_int = dbpath + pc_folder + "PathwayCommons.8.All.BINARY_SIF.hgnc.txt.sif";
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

    static String tp_ent = dbpath + transpath_folder + "reaction_semantic_full_entities_EDT.txt";
    static String tp_int = dbpath + transpath_folder + "reaction_semantic_full_connections_EDT.txt";
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

    static final String apop = "apoptosis" + sep;
    static final String inputAp = "input" + sep;
    static final String outputAp = "output" + sep;
    static final String ppiAp = "ppi" + sep;
    static final String outputFAp = syspath + apop + outputAp;
    static final String randomFAp = syspath + apop + outputAp + "random" + sep;
    static String hlAp = syspath + apop + inputAp + "HitGenes.txt";
    static String fpAp_dffAB = syspath + apop + inputAp + "FinalPlayers_dffAB.txt";
    static String fpAp_ppi = syspath + apop + inputAp + "FinalPlayers_ppi.txt";
    static String fpAp = syspath + apop + inputAp + "FinalPlayers.txt";
    static String foundDFAp = syspath + apop + outputAp + "foundDF_18_01.txt";
    static String foundDFAp_ppi = syspath + apop + outputAp + ppiAp + "foundDF.txt";
    static String foundDFAp_overrepr = syspath + apop + outputAp + "foundDF_overrepr.txt";
    static String foundDFAp_mirnas = syspath + apop + outputAp + "foundDFmirnas.txt";
    static String Aprandom_path = syspath + apop + outputAp + "random" + sep;
    static String Aprandom_path_ppi = syspath + apop + outputAp + ppiAp + "random" + sep;

    static final String apop_test = "apoptosis_test" + sep;
    static final String inputAp_test = "input" + sep;
    static final String outputAp_test = "output" + sep;
    static final String ppiAp_test = "ppi" + sep;
    static final String outputFAp_test = syspath + apop_test + outputAp;
    static final String randomFAp_test = syspath + apop_test + outputAp + "random" + sep;
    static String hlAp_test = syspath + apop_test + inputAp + "HitGenes.txt";
    static String fpAp_dffAB_test = syspath + apop_test + inputAp + "FinalPlayers_dffAB.txt";
    static String fpAp_ppi_test = syspath + apop_test + inputAp + "FinalPlayers_ppi.txt";
    static String fpAp_test = syspath + apop_test + inputAp + "FinalPlayers.txt";
    static String foundDFAp_test = syspath + apop_test + outputAp + "foundDF_18_01.txt";
    static String foundDFAp_ppi_test = syspath + apop_test + outputAp + ppiAp + "foundDF.txt";
    static String foundDFAp_overrepr_test = syspath + apop_test + outputAp + "foundDF_overrepr.txt";
    static String foundDFAp_mirnas_test = syspath + apop_test + outputAp + "foundDFmirnas.txt";
    static String Aprandom_path_test = syspath + apop_test + outputAp + "random" + sep;
    static String Aprandom_path_ppi_test = syspath + apop_test + outputAp + ppiAp + "random" + sep;

    static final String apop_test2 = "apoptosis_test2" + sep;
    static final String inputAp_test2 = "input" + sep;
    static final String outputAp_test2 = "output" + sep;
    static final String ppiAp_test2 = "ppi" + sep;
    static final String outputFAp_test2 = syspath + apop_test2 + outputAp;
    static final String randomFAp_test2 = syspath + apop_test2 + outputAp + "random" + sep;
    static String hlAp_test2 = syspath + apop_test2 + inputAp + "HitGenes.txt";
    static String fpAp_dffAB_test2 = syspath + apop_test2 + inputAp + "FinalPlayers_dffAB.txt";
    static String fpAp_ppi_test2 = syspath + apop_test2 + inputAp + "FinalPlayers_ppi.txt";
    static String fpAp_test2 = syspath + apop_test2 + inputAp + "FinalPlayers.txt";
    static String foundDFAp_test2 = syspath + apop_test2 + outputAp + "foundDF_18_01.txt";
    static String foundDFAp_ppi_test2 = syspath + apop_test2 + outputAp + ppiAp + "foundDF.txt";
    static String foundDFAp_overrepr_test2 = syspath + apop_test2 + outputAp + "foundDF_overrepr.txt";
    static String foundDFAp_mirnas_test2 = syspath + apop_test2 + outputAp + "foundDFmirnas.txt";
    static String Aprandom_path_test2 = syspath + apop_test2 + outputAp + "random" + sep;
    static String Aprandom_path_ppi_test2 = syspath + apop_test2 + outputAp + ppiAp + "random" + sep;

    static final String mapk = "MAPK" + sep;
    static final String inputMK = "input" + sep;
    static final String outputMK = "output" + sep;
    static final String ppiMK = "ppi" + sep;
    static final String outputFMK = syspath + mapk + outputMK;
    static final String randomFMK = syspath + mapk + outputMK + "random" + sep;
    static String hlMK = syspath + mapk + inputMK + "HitGenes.txt";
    static String fpMK_ppi = syspath + mapk + inputMK + "FinalPlayers_ppi.txt";
    static String fpMK = syspath + mapk + inputMK + "FinalPlayers.txt";
    static String foundDFMK = syspath + mapk + outputMK + "foundDF.txt";
    static String foundDFMK_ppi = syspath + mapk + outputMK + ppiMK + "foundDF.txt";
    static String foundDFMK_overrepr = syspath + mapk + outputMK + "foundDF_overrepr.txt";
    static String foundDFMK_mirnas = syspath + mapk + outputMK + "foundDFmirnas.txt";
    static String MKrandom_path = syspath + mapk + outputMK + "random" + sep;
    static String MKrandom_path_ppi = syspath + mapk + outputMK + ppiMK + "random" + sep;

    static final String mapk2 = "MAPK2" + sep;
    static final String inputMK2 = "input" + sep;
    static final String outputMK2 = "output" + sep;
    static final String ppiMK2 = "ppi" + sep;
    static final String outputFMK2 = syspath + mapk2 + outputMK2;
    static final String randomFMK2 = syspath + mapk2 + outputMK2 + "random" + sep;
    static String hlMK2 = syspath + mapk2 + inputMK + "HitGenes.txt";
    static String fpMK2_ppi = syspath + mapk2 + inputMK2 + "FinalPlayers_ppi.txt";
    static String fpMK2 = syspath + mapk2 + inputMK2 + "FinalPlayers.txt";
    static String foundDFMK2 = syspath + mapk2 + outputMK2 + "foundDF.txt";
    static String foundDFMK2_ppi = syspath + mapk2 + outputMK2 + ppiMK2 + "foundDF.txt";
    static String foundDFMK2_overrepr = syspath + mapk2 + outputMK2 + "foundDF_overrepr.txt";
    static String foundDFMK2_mirnas = syspath + mapk2 + outputMK2 + "foundDFmirnas.txt";
    static String MK2random_path = syspath + mapk2 + outputMK2 + "random" + sep;
    static String MK2random_path_ppi = syspath + mapk2 + outputMK2 + ppiMK2 + "random" + sep;

    static final String pathways = "pathways" + sep;
    static final String inputPTH = "input" + sep;
    static final String outputPTH = "output" + sep;
    static final String outputFPTH = syspath + pathways + outputPTH;
    static final String randomFPTH = syspath + pathways + outputPTH + "random" + sep;
    static String hlPTH = syspath + pathways + inputPTH + "HitGenes.txt";
    static String fpPTH = syspath + pathways + inputPTH + "FinalPlayers.txt";
    static String foundDFPTH = syspath + pathways + outputPTH + "foundDF.txt";
    static String folder_foundDFPTH_hubs = syspath + pathways + outputPTH + "hubs" + sep;
    static String foundDFPTH_overrepr = syspath + pathways + outputPTH + "foundDF_overrepr.txt";
    static String PTHrandom_path = syspath + pathways + outputPTH + "random" + sep;

    static final String adnc = "Adenocarcinoma" + sep;
    static final String inputAd = "input" + sep;
    static final String outputAd = "output" + sep;
    static final String outputFAd = syspath + adnc + outputAd;
    static final String outputFAd_ntwrk_dir = syspath + adnc + outputAd + "networks" + sep;
    static String hlAd = syspath + adnc + inputAd + "HitGenes.txt";
    static String fpAd = syspath + adnc + inputAd + "FinalPlayers.txt";
    static String foundDFAd = syspath + adnc + outputAd + "foundDF.txt";
    static String foundDFAd_overrepr = syspath + adnc + outputAd + "foundDF_overrepr.txt";
    static String foundDFAd_mirnas = syspath + adnc + outputAd + "foundDFmirnas.txt";
    static String Adrandom_path = syspath + adnc + outputAd + "random" + sep;
    static final String ppiAd = "ppi" + sep;
    static String foundDFAd_ppi = syspath + adnc + outputAd + ppiAd + "foundDF.txt";
    static String Adrandom_path_ppi = syspath + adnc + outputAd + ppiAd + "random" + sep;
    static String fpAd_ppi = syspath + adnc + inputAd + "FinalPlayers_ppi.txt";

    static final String arp = "ARP" + sep;
    static final String inputARP = "input" + sep;
    static final String outputARP = "output" + sep;
    static final String outputFARP = syspath + arp + outputARP;
    static final String outputFARP_ntwrk_dir = syspath + arp + outputARP + "networks" + sep;
    static String hlARP = syspath + arp + inputARP + "HitGenes.txt";
    static String hlARP_LIMCH1 = syspath + arp + inputARP + "HitGenes_LIMCH1.txt";
    static String fpARP = syspath + arp + inputARP + "FinalPlayers.txt";
    static String foundDFARP = syspath + arp + outputARP + "foundDF.txt";
    static String foundDFARP_overrepr = syspath + arp + outputARP + "foundDF_overrepr.txt";
    static String foundDFARP_mirnas = syspath + arp + outputARP + "foundDFmirnas.txt";
    static String ARPrandom_path = syspath + arp + outputARP + "random" + sep;
    static final String ppiARP = "ppi" + sep;
    static String foundDFARP_ppi = syspath + arp + outputARP + ppiARP + "foundDF.txt";
    static String ARPrandom_path_ppi = syspath + arp + outputARP + ppiARP + "random" + sep;
    static String fpARP_ppi = syspath + arp + inputARP + "FinalPlayers_ppi.txt";

    static final String muscles = "Muscle Differentiation" + sep;
    static final String input = "input" + sep;
    static final String output = "output" + sep;
    static final String outputF = syspath + muscles + output;
    static String hl = syspath + muscles + input + "HitGenes.txt";
    static String fp = syspath + muscles + input + "FinalPlayers.txt";
    static String foundDF = syspath + muscles + output + "foundDF_28_11.txt";
    static String foundDF_overrepr = syspath + muscles + output + "foundDF_overrepr.txt";
    static String foundDF_mirnas = syspath + muscles + output + "foundDFmirnas.txt";

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

    static final String case2e = "Case2_Control_Early" + sep;
    static final String inputC2e = "input" + sep;
    static final String outputC2e = "output" + sep;
    static final String outputFC2E = syspath + case2e + outputC2e;
    static String hlC2e = syspath + case2e + inputC2e + "HitGenes.txt";
    static String fpC2e = syspath + case2e + inputC2e + "FinalPlayers.txt";
    static String exprC2e = syspath + case2e + inputC2e + "expression.txt";
    static String foundDFC2e = syspath + case2e + outputC2e + "foundDF.txt";
    static String foundDFC2e_overrepr = syspath + case2e + outputC2e + "foundDF_overrepr.txt";
    static String foundDFC2e_mirnas = syspath + case2e + outputC2e + "foundDFmirnas.txt";

    static final String case3 = "Case3_Hypertrophy" + sep;
    static final String inputC3 = "input" + sep;
    static final String outputC3 = "output" + sep;
    static final String outputFC3 = syspath + case3 + outputC3;
    static String hlC3 = syspath + case3 + inputC3 + "HitGenes.txt";
    static String fpC3 = syspath + case3 + inputC3 + "FinalPlayers.txt";
    static String foundDFC3 = syspath + case3 + outputC3 + "foundDF_18_01.txt";
    static String foundDFC3_overrepr = syspath + case3 + outputC3 + "foundDF_overrepr.txt";
    static String foundDFC3_mirnas = syspath + case3 + outputC3 + "foundDFmirnas.txt";

    static String muscle_nw2 = syspath + muscles + output + "nw2_edt.txt";
    static String muscle_nw2_w_connectivities = syspath + muscles + output + "nw2_edt_con.txt";

    static final String stat = "Statistics" + sep;
    static final String inputstatf = "input" + sep;
    static final String outputstatf = "output" + sep;
    static final String outputstat = syspath + stat + outputstatf;

}
