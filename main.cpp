#include <bits/stdc++.h>

using namespace std;

int sizeseq(char seq[]){
    string s=seq;
    return s.length();
}

class Sequence
{
  	protected:
        char * seq;
    public:
 	 	// constructors and destructor
 	 	// pure virtual function that should be overridden because every
        // type of sequence has its own details to be printed
        virtual void Print()= 0;
 	 	// friend function that will find the LCS (longest common
        // subsequence) between 2 sequences of any type, according to
        // polymorphism
        friend char* Align(Sequence * s1, Sequence * s2);
};
enum DNA_Type{promoter, motif, tail, noncoding};
class DNA : public Sequence
{
  	private:
        DNA_Type type;
        DNA * complementary_strand;
        int startIndex;
        int endIndex;
    public:
 	 	// constructors and destructor
        DNA(){}
        DNA(char * seqq, DNA_Type atype){
            seq=new char[sizeseq(seqq)];
            for(int i=0;i<sizeseq(seqq);i++){
                seq[i]=seqq[i];
            }
            type=atype;
        }
        DNA(const DNA& rhs){
            seq=new char[sizeseq(rhs.seq)];
            for(int i=0;i<sizeseq(rhs.seq);i++){
                seq[i]=rhs.seq[i];
            }
            type=rhs.type;
        }
 	 	// function printing DNA sequence information to user
        void Print(){
            cout<<seq;
        }
        // function to convert the DNA sequence to RNA sequence
        // It starts by building the complementary_strand of the current
        // DNA sequence (starting from the startIndex to the endIndex), then,
        // it builds the RNA corresponding to that complementary_strand.
//        RNA ConvertToRNA(){
//            char* b;
//            BuildComplementaryStrand();
//            for(int i=0;i<endIndex-startIndex;i++){
//                if(complementary_strand->seq[i]=='T')b[i]='U';
//                else
//                    b[i]=complementary_strand->seq[i];
//                if(i==endIndex-startIndex-1)b[i+1]=NULL;
//            }
//            return RNA(b,mRNA);
//        }
        // function to build the second strand/pair of DNA sequence
	    // To build a complementary_strand (starting from the startIndex to
        // the endIndex), convert each A to T, each T to A, each C to G, and
        // each G to C. Then reverse the resulting sequence.
        void BuildComplementaryStrand(){
            for(int i=startIndex;i<endIndex;i++){
                if(seq[i]=='A')complementary_strand->seq[i]='T';
                else if(seq[i]=='T')complementary_strand->seq[i]='A';
                else if(seq[i]=='C')complementary_strand->seq[i]='G';
                else if(seq[i]=='G')complementary_strand->seq[i]='C';
            }
            reverse(complementary_strand->seq,complementary_strand->seq+endIndex);
        }
  };

// struct representing a codon of 3 DNA/RNA characters and ‘\0’
struct Codon
{
  	char value[4];    	    // 4th location for null character
    char AminoAcid;  	    // corresponding AminoAcid according to Table
};

// need to create one object of that class to load the AminoAcids table
// into memory
class CodonsTable
{
  	private:
        Codon codons[64];
   	public:
 	 	// function to load all codons from the given text file
        void LoadCodonsFromFile(){
        ifstream thefile("codonsfile.txt");
        char value[4];
        char AminoAcid;
        int i=0;
        while(thefile>>value>>AminoAcid){
            setCodon(value,AminoAcid,i);
            i++;
        }
        thefile.close();

        }
        Codon getAminoAcid(char value[4]){
            int test=0;
            for(int i=0;i<64;i++){
                for(int j=0;j<4;j++){
                    if(codons[i].value[j]==value[j])test++;
                    if (test==4)return codons[i];
                }
                test=0;

            }

        }
        void setCodon(char value[4], char AminoAcid, int index){
            for(int i=0;i<4;i++){
                    codons[index].value[i]=value[i];
                }
            codons[index].AminoAcid=AminoAcid;
        }
        void gettest(){
        for(int i=0;i<64;i++){cout<<codons[i].value<<' '<<codons[i].AminoAcid<<endl;}

        }
};

enum RNA_Type {mRNA, pre_mRNA, mRNA_exon, mRNA_intron};
class RNA : public Sequence
{
  	private:
        RNA_Type type;
  	public:
 	 	// constructors and destructor
        RNA();
        RNA(char * seqq, RNA_Type atype){
            seq=new char[sizeseq(seqq)];
            for(int i=0;i<sizeseq(seqq);i++){
                seq[i]=seqq[i];
            }
            type=atype;
        }
        RNA(const RNA& rhs){
            seq=new char[sizeseq(rhs.seq)];
            for(int i=0;i<sizeseq(rhs.seq);i++){
                seq[i]=rhs.seq[i];
            }
            type=rhs.type;
        }
 	 	// function to be overridden to print all the RNA information
        void Print(){
            cout<<seq<<endl;
        }
 	 	// function to convert the RNA sequence into protein sequence
        // using the codonsTable object
//        Protein ConvertToProtein(CodonsTable & table){
//            char* a=new char[3];
//            table.LoadCodonsFromFile();
//            char* b=new char[sizeseq(seq)/3];
//            for(int i=0;i<sizeseq(seq);i+=3){
//                a[0]=seq[i];
//                a[1]=seq[i+1];
//                a[2]=seq[i+2];
//                if(table.getAminoAcid(a).AminoAcid!='0')
//                    b[i/3]=table.getAminoAcid(a).AminoAcid;
//            }
//        }
 	 	// function to convert the RNA sequence back to DNA
        DNA ConvertToDNA(){
            char* b;
            int sizeA=sizeseq(seq);
            b=new char[sizeA];
            for(int i=0;i<=sizeA;i++){
                if(seq[i]=='U')b[i]='A';
                else if(seq[i]=='A')b[i]='T';
                else if(seq[i]=='G')b[i]='C';
                else if(seq[i]=='C')b[i]='G';
            }
            reverse(b,b+sizeA);
            DNA a(b,promoter);
            cout<<"";
            return a;
        }
};

enum Protein_Type {Hormon, Enzyme, TF, Cellular_Function};

class Protein : public Sequence
{
  	private:
        Protein_Type type;
      public:
 	 	// constructors and destructor
 	 	Protein();
 	 	Protein(char * p){
            seq=new char[sizeseq(p)];
            for(int i=0;i<sizeseq(p);i++){
                seq[i]=p[i];
            }
 	 	}
 	 	// return an array of DNA sequences that can possibly
        // generate that protein sequence
        DNA GetDNAStrandsEncodingMe(DNA bigDNA){
        }
};
int main()
{
    int n;
    cin>>n;
    char* arr=new char[n];
    cin>>arr;
    RNA a(arr,mRNA);
    DNA b=a.ConvertToDNA();
    b.Print();
    return 0;
}

