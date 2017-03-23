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

static void generate_input (const Dataset &data);
static void run_R (const Dataset &data);

static void generateRplot(Dataset &data){

    generate_input(data);
    run_R(data);
}

static  void generate_input(const Dataset &data){

    stringstream Rout;
    Rout << data.output_dir  << "/generateRinput.input";
    ofstream Rfile;
    Rfile.open(Rout.str());

    //Not sure if the original print function prints with spaces
    //between commands or newline characters so for now I have newline
    //characters.
    // REPLY: I believe this is correct, If I understand the "<<" command or option
    // or whatever it is correctly, then the text is printed literally as it appears
    Rfile << "pdf(\"" << data.output_dir << "" << data.TF_name << "_semplot.pdf\")\n"
          <<"source(\"src/plotSEM_Functions.R\")\n"
          <<"plotSEM(\"" << data.output_dir << "\", \"" << data.TF_name << "\", error=TRUE)\n"
          <<"dev.off ()";
    Rfile.close();
}

static void run_R(const Dataset &data){

    stringstream Rout;
    Rout << "R --vanilla < " << data.output_dir << "/generateRinput.input";
    system(Rout.str().c_str());
    Rout.clear();
    Rout << "rm " << data.output_dir << "/generateRinput.input";
    system(Rout.str().c_str());

}
