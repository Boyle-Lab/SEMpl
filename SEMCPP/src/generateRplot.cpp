#include "iterativeSEM.hpp"
#include <sstream>
#include <stdlib.h>
#include <fstream>

/*
*
*   Effect: Creates the input for the R plot and runs the
*   plotSEM_Functions file.
*
*/

using namespace std;

void generate_input (Dataset &data);
void run_R (Dataset &data);

void generateRplot(Dataset &data){

    generate_input(data);
    run_R(data);
}

void generate_input(Dataset &data){

    stringstream Rout;
    Rout << data.output_dir  << "/generateRinput.input";
    ofstream Rfile;
    Rfile.open(Rout.str());

    //Not sure if the original print function prints with spaces
    //between commands or newline characters so for now I have newline
    //characters.
    Rfile << "pdf(\"$output/$TFname_semplot.pdf\")\n"
          <<"source(\"src/plotSEM_Functions.R\")\n"
          <<"plotSEM(\"$output\", \"$TFname\", error=TRUE)\n"
          <<"dev.off ()";
    Rfile.close();
}

void run_R(Dataset &data){

    stringstream Rout;
    Rout << "R --vanilla < " << data.output_dir << "/generateRinput.input";
    system(Rout.str().c_str());
    Rout.clear();
    Rout << "rm " << data.output_dir << "/generateRinput.input";
    system(Rout.str().c_str());

}
