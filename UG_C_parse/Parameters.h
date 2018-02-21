//
// Created by swoo on 9/15/2017.
//
//#ifndef UGRAPHS_PARAMETERS_H
//#define UGRAPHS_PARAMETERS_H
//
//#endif //UGRAPHS_PARAMETERS_H
#include <string>

class Param{
    private:
        static Param* instance;
        int MAX_IDX = 110;
        int MIN_IDX = 10;
        string fasta_file = "H:\\Ensembl\\unmasked\\Homo_sapiens.GRCh38.dna.chromosome.1.strip.fa";
        string test_file = "H:\\CLionProjects\\UGraphs\\Homo_sapiens.GRCh38.dna.chromosome.1.test.txt";
        string gff_file = "H:\\Ensembl\\GFF\\Homo_sapiens.GRCh38.88.chr1.gff3";
//        string gff_file = "H:\\Ensembl\\GFF\\test1.txt";
//        string test_file = "H:\\Ensembl\\GFF\\test_out.txt";
        string chr_in_process = "1";

        Param(); /* Private constructor to prevent instancing. */

    public:
        static Param* getInstance(); /* Static access method. */
        int get_MAX_IDX(){return MAX_IDX;}
        int get_MIN_IDX(){return MIN_IDX;}
        string get_fasta_file(){return fasta_file;}
        string get_test_file(){return test_file;}
        string get_gff_file(){return gff_file;}
        string get_chr_in_process(){return chr_in_process;}
};

Param* Param::instance = 0; /* Null, because instance will be initialized on demand. */
Param* Param::getInstance(){
    if (instance == 0){instance = new Param();}
    return instance;
}

Param::Param()
{}
