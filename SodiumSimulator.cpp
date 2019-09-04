#include <bits/stdc++.h>

using namespace std;


#define SCALE 0.001
#define INV_SCALE 1/SCALE

#define DIM_X 5
#define DIM_Y 5
#define DIM_Z 3
#define ALPHA 0.01
#define PI 3.14159265359

double conc, tx_conc, diameter_cell = 5, deltaamp;

// #### GAP JUNCTIONS PROBABILITIES FOR ASTROCYTES #### //
double phl[50] = {3.33333333333e-01,9.51756626745e-02,2.71812917035e-02,7.76661124715e-03,2.22288659849e-03,6.39921679991e-04,1.87410416804e-04,5.81750667195e-05,2.17596143238e-05,1.1436923158e-05,7.88454209682e-06,7.43738619183e-06,7.37970057786e-06,7.29603316347e-06,7.27478942971e-06,7.26006289992e-06,7.26084787208e-06,7.26080132601e-06,7.26061054996e-06,7.26081361742e-06,7.26079620991e-06,7.26072365567e-06,7.26058079345e-06,7.26074725419e-06,7.26087576894e-06,7.26073008288e-06,7.26061194028e-06,7.26074336727e-06,7.26101528686e-06,7.26081974085e-06,7.26091667847e-06,7.26058059059e-06,7.26084014577e-06,7.2610063969e-06,7.26069065682e-06,7.26083092741e-06,7.26076153595e-06,7.26071756287e-06,7.26092535023e-06,7.26076421324e-06,7.26060026219e-06,7.26075209967e-06,7.26093367537e-06,7.26073986493e-06,7.26039032094e-06,7.26091299989e-06,7.26077756319e-06,7.26071491915e-06,7.2607710224e-06,7.26082337127e-06};
double plh[50] = {0.333333333333,0.706523083189,0.825908352326,0.86117728508,0.871356970354,0.874272895353,0.875107602476,0.875346523333,0.875414979813,0.875433066086,0.875439533589,0.875437420808,0.875439605921,0.875443525702,0.875437685437,0.875440602592,0.875441218644,0.875440938251,0.875438148928,0.875441400815,0.875441147657,0.875440025483,0.875437801975,0.875440355641,0.875442338151,0.875440081279,0.87543825197,0.875440285023,0.87544449347,0.875441466845,0.875442967173,0.875437765288,0.875441782611,0.875444355802,0.875439468868,0.875441639938,0.875440565916,0.875439885313,0.875443101387,0.875440607355,0.875438069767,0.875440419864,0.875443230241,0.875440230498,0.875434820356,0.875442910232,0.875440813981,0.875439844395,0.875440712745,0.875441522986};
double phh[50] = {0.333333333333,0.198301254137,0.14691035597,0.131056103673,0.126420143048,0.125087182967,0.124704987107,0.1245953016,0.124563260573,0.124555496991,0.124552581869,0.124555141806,0.124553014379,0.124549178265,0.124555039774,0.124552137345,0.124551520508,0.124551800947,0.124554590462,0.124551338371,0.124551591547,0.124552713793,0.124554937445,0.124552383612,0.124550400973,0.124552657991,0.124554487418,0.124552454233,0.124548245514,0.124551272335,0.12454977191,0.124554974131,0.124550956549,0.124548383192,0.124553270441,0.124551099232,0.124552173322,0.124552853969,0.124549637687,0.124552131881,0.124554669633,0.124552319384,0.124549508825,0.124552508762,0.124557919254,0.124549828855,0.124551925241,0.12455289489,0.124552026484,0.124551216191};

// Classe para representar todas as variáveis da célula (Ex.: astrócito)
class Cell {
public:
	map<string, double> parameters;
	int id;

	Cell() {
		/*
			0v1: taxa de liberação de Cálcio do Tx
			1Y:
			2vin: o fluxo de calcio que parte do espaço extracelular, por meio da membrana do astrócito, até o interior do citosol
			3VM2: o fluxo máximo de íons de cálcio fora da bomba (???)
			4C: Concentração de Cálcio no citosol
			5n: Coeficiente de Hill (2.02)
			6K2:
			10kf: constante que determina a liberação de cálcio do RE para o citosol
			22D: coneficiente de difusão
			23l: volume da célula
			24K: taxa máxima da ativação do receptor (Nakano, 2010; Eq. 3)
			25ka: taxa máxima da ativação do receptor (2.5) (Nakano, 2010; Eq. 3)
			26m: Coeficiente de Hill (2.2)
			27phh - 29plh: probabilidades das gap junctions
		*/

		double v1 = 100*SCALE; parameters["v1"] = v1;
		double Y = 0*SCALE; parameters["Y"] = Y;
		double vin = 0.05*SCALE; parameters["vin"] = vin;
		double VM2 = 15*SCALE; parameters["VM2"] = VM2;
		double C = 0.1*SCALE; parameters["C"] = C;
		double n = 2.02*SCALE; parameters["n"] = n;
		double K2 = 0.1*SCALE; parameters["K2"] = K2;
		double VM3 = 40*SCALE; parameters["VM3"] = VM3; // Porque nao 40, como diz no artigo?
		double ko = 0.5*SCALE; parameters["ko"] = ko;
		double ER = 1.5*SCALE; parameters["ER"] = ER;
		double kf = 0.5*SCALE; parameters["kf"] = kf;
		double kp = 0.3*SCALE; parameters["kp"] = kp;
		double kdeg = 0.08*SCALE; parameters["kdeg"] = kdeg;
		double vp = 0.05*SCALE; parameters["vp"] = vp;
		double kcaaa = 0.15*SCALE; parameters["kcaaa"] = kcaaa;
		double kcai = 0.15*SCALE; parameters["kcai"] = kcai;
		double kip3 = 0.1*SCALE; parameters["kip3"] = kip3;
		double IP3 = 0.1*SCALE; parameters["IP3"] = IP3;
		double q = 4*SCALE; parameters["q"] = q;
		double W = 0*SCALE; parameters["W"] = W;
		double A = 0*SCALE; parameters["A"] = A;
		double kia = 0.5*SCALE; parameters["kia"] = kia;
		double D = 350 * (4 * PI * diameter_cell / 3.0)/*122500*/; parameters["D"] = D;
		double l = PI *SCALE * pow(diameter_cell / 2, 2); parameters["l"] = l;
		double K = 0.0006 *SCALE; parameters["K"] = K;
		double ka = 2.5*SCALE; parameters["ka"] = ka;
		double m = 2.2*SCALE; parameters["m"] = m;
		parameters["phl"] = phl[0];
		parameters["plh"] = plh[0];
		parameters["phh"] = phh[0];

		double Na_i = SCALE*15000;  parameters["Na_i"] = Na_i; //Langer3, Chatton 2016 - Colocar uma relação
		double Na_o = SCALE*150000/(DIM_X * DIM_Y * DIM_Z);  parameters["Na_o"] = Na_o; //chatton 2016
		double NaD = (4/3)*PI*diameter_cell*6000; parameters["NaD"] = NaD; 
		double C_o = SCALE*2300/(DIM_X * DIM_Y * DIM_Z); parameters["C_o"] = C_o; //Kirischuk 1997
		double Vm = -80; parameters["Vm"] = Vm; //Verify this value later - Kirischuk 2012

	}

	void setId(int ID) {
		id = ID;
		parameters["id"] = ID;
	}
};

class Network {
public:
	int NC; // Number of cells
	Cell tecido[DIM_Y][DIM_X][DIM_Z];
	list<int> *connect;

	Network() {
		NC = DIM_X * DIM_Y * DIM_Z;
		connect = new list<int>[NC];

		int id_cont = 0;

		for (int k = 0; k < DIM_Z; k++) {
			for (int i = 0; i < DIM_Y; i++) {
				for (int j = 0; j < DIM_X; j++) {
					Cell celula;
					celula.setId(id_cont);
					tecido[i][j][k] = celula;
					id_cont++;
				}
			}
		}
	}

	int getId(int x, int y, int z) {
		return tecido[y][x][z].parameters["id"];
	}

	double get(int x, int y, int z, string parameter) {
		return tecido[y][x][z].parameters[parameter];
	}

	double get(int id, string parameter) {
		int x = id % DIM_X;
		int y = (id / DIM_X) % DIM_Y;
		int z = id / (DIM_X * DIM_Y);
		return tecido[y][x][z].parameters[parameter];
	}

	void set(int x, int y, int z, string parameter, double value) {
		tecido[y][x][z].parameters[parameter] = value;
	}

	void set(int id, string parameter, double value) {
		int x = id % DIM_X;
		int y = (id / DIM_X) % DIM_Y;
		int z = id / (DIM_X * DIM_Y);
		tecido[y][x][z].parameters[parameter] = value;
	}

	void accumulate(int x, int y, int z, string parameter, double add) {
		tecido[y][x][z].parameters[parameter] += add;
	}

	void accumulate(int id, string parameter, double add) {
		int x = id % DIM_X;
		int y = (id / DIM_X) % DIM_Y;
		int z = id / (DIM_X * DIM_Y);
		tecido[y][x][z].parameters[parameter] += add;
	}

	void changeSignal(int x, int y, int z, string parameter) {
		tecido[y][x][z].parameters[parameter] *= -1;
	}

	void printTissue() {
		cout << endl;
		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				cout << fixed << setprecision(2) << get(i, j, trunc(DIM_Z / 2), "Na_i") << "  ";
			}
			cout << endl;
		}
	}

	void writeFileHeader(ofstream &file){
		for (int i = 0; i < DIM_Y; i++) {
		  for (int j = 0; j < DIM_X; j++) {
			if(i == DIM_Y-1 && j == DIM_X-1)
			  file << i << "_" << j;
			else
			  file <<  i << "_" << j << ",";
		  }
		}
		file << endl;
	}

	void printTissuef(ofstream &file, int time_int){
		file << time_int << ",";
		for (int i = 0; i < DIM_Y; i++) {
		  for (int j = 0; j < DIM_X; j++) {
			if(i == DIM_Y-1 && j == DIM_X-1)
			  file << get(j, i, trunc(DIM_Z / 2), "C");
			else
			  file <<  get(j, i, trunc(DIM_Z / 2), "C") << ",";
		  }
		}
		file << endl;
	}

	void regularDegree() {
		for (int x = 0; x < DIM_X; x++) {
			for (int y = 0; y < DIM_Y; y++) {
				for (int z = 0; z < DIM_Z; z++) {
					// x + 1
					if (y + 1 < DIM_Y) {
						connect[getId(x, y, z)].push_back(getId(x, y + 1, z));
					} else {
						connect[getId(x, y, z)].push_back(-1);
					}
					// x - 1
					if (y - 1 >= 0) {
						connect[getId(x, y, z)].push_back(getId(x, y - 1, z));
					} else {
						connect[getId(x, y, z)].push_back(-1);
					}
					// y + 1
					if (x + 1 < DIM_X) {
						connect[getId(x, y, z)].push_back(getId(x + 1, y, z));
					} else {
						connect[getId(x, y, z)].push_back(-1);
					}
					// y - 1
					if (x - 1 >= 0) {
						connect[getId(x, y, z)].push_back(getId(x - 1, y, z));
					} else {
						connect[getId(x, y, z)].push_back(-1);
					}
					// z + 1
					if (z + 1 < DIM_Z) {
						connect[getId(x, y, z)].push_back(getId(x, y, z + 1));
					} else {
						connect[getId(x, y, z)].push_back(-1);
					}
					// z - 1
					if (z - 1 >= 0) {
						connect[getId(x, y, z)].push_back(getId(x, y, z - 1));
					} else {
						connect[getId(x, y, z)].push_back(-1);
					}
				}
			}
		}

		// IMPRIMINDO AS CONEXÕES DE CADA CÉLULA
		list<int>::iterator it;

		for (int v = 0; v < NC; v++) {
			cout << "Célula " << v << ": ";
			for (it = connect[v].begin(); it != connect[v].end(); it++){
				cout << *it << " ";
			}
			cout << endl;
		}
	}

	void linkRadius(int radius) {
		for (int x = 0; x < DIM_X; x++) {
			for (int y = 0; y < DIM_Y; y++) {
				for (int z = 0; z < DIM_Z; z++) {
					for (int r = 1; r <= radius; r++) {
						// x + r
						if (y + r < DIM_Y) {
							connect[getId(x, y, z)].push_back(getId(x, y + r, z));
						} else {
							connect[getId(x, y, z)].push_back(-1);
						}
						// x - r
						if (y - r >= 0) {
							connect[getId(x, y, z)].push_back(getId(x, y - r, z));
						} else {
							connect[getId(x, y, z)].push_back(-1);
						}
						// y + r
						if (x + r < DIM_X) {
							connect[getId(x, y, z)].push_back(getId(x + r, y, z));
						} else {
							connect[getId(x, y, z)].push_back(-1);
						}
						// y - r
						if (x - r >= 0) {
							connect[getId(x, y, z)].push_back(getId(x - r, y, z));
						} else {
							connect[getId(x, y, z)].push_back(-1);
						}
						// z + r
						if (z + r < DIM_Z) {
							connect[getId(x, y, z)].push_back(getId(x, y, z + r));
						} else {
							connect[getId(x, y, z)].push_back(-1);
						}
						// z - r
						if (z - r >= 0) {
							connect[getId(x, y, z)].push_back(getId(x, y, z - r));
						} else {
							connect[getId(x, y, z)].push_back(-1);
						}
					}
				}
			}
		}

		// IMPRIMINDO AS CONEXÕES DE CADA CÉLULA
		list<int>::iterator it;

		for (int v = 0; v < NC; v++) {
			cout << "Célula " << v << ": ";
			for (it = connect[v].begin(); it != connect[v].end(); it++){
				cout << *it << " ";
			}
			cout << endl;
		}
	}

	vector<int> getConnections(int id) {
		vector<int> connections;
		list<int>::iterator it;

		for (it = connect[id].begin(); it != connect[id].end(); it++){
			connections.push_back(*it);
		}

		return connections;
	}

	int numberConnections() {
		return connect[0].size();
	}
};

class Gillespie {
	Network *tecido;

public:
	Gillespie(Network *rede) {
		tecido = rede;

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				cout << "(" << j << "," << i << ") " << tecido->getId(i, j, trunc(DIM_Z / 2)) << " " << tecido->get(i, j, trunc(DIM_Z / 2), "C") << " ";
			}
			cout << endl;
		}

		//cout << sigma0(0) << endl;
	}

	// Reaction 1
	double sigma0(int id) {
		return tecido->get(id, "vin");
	}

	// Reaction 2
	double sigma1(int id) {
		return 4 * tecido->get(id, "VM3") * ((pow(tecido->get(id, "kcaaa"), tecido->get(id, "n")) * pow(tecido->get(id, "C"), tecido->get(id, "n")) ) / ((pow(tecido->get(id, "C"), tecido->get(id, "n")) + pow(tecido->get(id, "kcaaa"), tecido->get(id, "n"))) * (pow(tecido->get(id, "C"), tecido->get(id, "n")) + pow(tecido->get(id, "kcai"), tecido->get(id, "n"))) ) * (pow(tecido->get(id, "IP3"), tecido->get(id, "m")) / (pow(tecido->get(id, "kip3"), tecido->get(id, "m")) + pow(tecido->get(id, "IP3"), tecido->get(id, "m")) ) )) * (tecido->get(id, "ER") - tecido->get(id, "C"));
	}

	// Reaction 3
	double sigma2(int id) {
		return tecido->get(id, "VM2") * ( pow(tecido->get(id, "C"), 2) / (pow(tecido->get(id, "K2"), 2) + pow(tecido->get(id, "C"), 2)) );
	}

	// Reaction 4
	double kf_Ea(int id) {
		return tecido->get(id, "kf") * tecido->get(id, "ER");
	}

	// Reaction 5
	double kf_Ca(int id) {
		return tecido->get(id, "kf") * tecido->get(id, "C");
	}

	// Reaction 6
	double ko_Ca(int id) {
		return tecido->get(id, "ko") * tecido->get(id, "C");
	}

	// Reaction 7
	double sigma3(int id) {
		return tecido->get(id, "vp") * (pow(tecido->get(id, "C"), 2) / (pow(tecido->get(id, "C"), 2) + pow(tecido->get(id, "kp"), 2)) );
	}

	// Reaction 8
	double kd_Ia(int id) {
		return tecido->get(id, "kdeg") * tecido->get(id, "IP3");
	}

	// Reaction 9
	vector<double> diffusions(int id) {
		int nConnections = tecido->numberConnections();
		vector<double> diffusions(nConnections * 3);

		vector<int> connections(nConnections);
		connections = tecido->getConnections(id);
		double value;

		for (int i = 0; i < connections.size(); i++){
			for (int gj = 0; gj < 3; gj++) {
				if (connections[i] != -1)
					value = diffusionEquation(id, connections[i], gj);
				else
					value = 0;

				diffusions[(3 * i) + gj] = value;
			}
		}

		return diffusions;
	}

	double diffusionEquation(int id1, int id2, int gap_junction) {
		double vol_cell = (4 / 3) * (PI * pow((diameter_cell / 2), 3));
		double diff;

		if (tecido->get(id1, "C") <= tecido->get(id2, "C"))
			return 0;

		if (gap_junction == 0) {
			diff = (tecido->get(id1, "D") / vol_cell) * (tecido->get(id1, "C") - tecido->get(id2, "C")) * tecido->get(id1, "phh");
		} else if (gap_junction == 1) {
			diff = (tecido->get(id1, "D") / vol_cell) * (tecido->get(id1, "C") - tecido->get(id2, "C")) * tecido->get(id1, "phl");
		} else {
			diff = (tecido->get(id1, "D") / vol_cell) * (tecido->get(id1, "C") - tecido->get(id2, "C")) * tecido->get(id1, "plh");
		}

		return diff;
	}

	// Reaction 10
	vector<double> Na_i_diffusions(int id) {
		int nConnections = tecido->numberConnections();
		vector<double> diffusions(nConnections * 3);

		vector<int> connections(nConnections);
		connections = tecido->getConnections(id);
		double value;

		for (int i = 0; i < connections.size(); i++){
			for (int gj = 0; gj < 3; gj++) {
				if (connections[i] != -1)
					value = Na_i_diffusionEquation(id, connections[i], gj);
				else
					value = 0;

				diffusions[(3 * i) + gj] = value;
			}
		}

		return diffusions;
	}

	double Na_i_diffusionEquation(int id1, int id2, int gap_junction) {
		double vol_cell = (4 / 3) * (PI * pow((diameter_cell / 2), 3));
		double diff;

		if (tecido->get(id1, "Na_i") <= tecido->get(id2, "Na_i"))
			return 0;

		if (gap_junction == 0) {
			diff = (tecido->get(id1, "NaD") / vol_cell) * (tecido->get(id1, "Na_i") - tecido->get(id2, "Na_i")) * tecido->get(id1, "phh");
		} else if (gap_junction == 1) {
			diff = (tecido->get(id1, "NaD") / vol_cell) * (tecido->get(id1, "Na_i") - tecido->get(id2, "Na_i")) * tecido->get(id1, "phl");
		} else {
			diff = (tecido->get(id1, "NaD") / vol_cell) * (tecido->get(id1, "Na_i") - tecido->get(id2, "Na_i")) * tecido->get(id1, "plh");
		}

		return diff;
	}

	vector<double> run() {
		int nConnections = tecido->numberConnections();
		int NC = DIM_X * DIM_Y * DIM_Z; // Total number of cells
		int num_reactions = 10; // Number of reactions (7 Intracellular + 2 Intercelular)  	/* Por que uma a mais? */
		double max_reaction = 0, reaction_choice, alfa_0 = 0, reaction_value;
		vector<double> retorno(5);
		vector<double> diffusion(nConnections * 3);
		vector<double> Na_i_diffusion(nConnections * 3);
		double reactions[DIM_Y][DIM_X][DIM_Z][num_reactions + 2*(nConnections * 3 - 1)]; //0 a 43

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {

				for (int k = 0; k < DIM_Z; k++) {
					// << Begin Reactions
					for (int r = 0; r < num_reactions; r++) {
						if (r == 0) {
							reaction_value = sigma0(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 1) {
							reaction_value = sigma1(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 2) {
							reaction_value = sigma2(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 3) {
							reaction_value = kf_Ea(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 4) {
							reaction_value = kf_Ca(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 5) {
							reaction_value = ko_Ca(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 6) {
							reaction_value = sigma3(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 7) {
							reaction_value = kd_Ia(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 8) {
							diffusion = diffusions(tecido->getId(i, j, k));
							for (int d = 0; d < diffusion.size(); d++) {
								reaction_value = diffusion[d];
								reactions[j][i][k][r + d] = reaction_value;

								/*if (reaction_value >= max_reaction) {
									max_reaction = reaction_value;
									cell[0] = j; cell[1] = i; cell[2] = k;

									if (d >= 0 && d <= 2)
										reaction_choice = 8;
									else if (d >= 3 && d <= 5)
										reaction_choice = 9;
									else if (d >= 6 && d <= 8)
										reaction_choice = 10;
									else if (d >= 9 && d <= 11)
										reaction_choice = 11;
									else if (d >= 12 && d <= 14)
										reaction_choice = 12;
									else
										reaction_choice = 13;
								}*/

								alfa_0 += reaction_value;
							}
						} else if (r == 9) {
							Na_i_diffusion = Na_i_diffusions(tecido->getId(i, j, k));
							for (int d = diffusion.size()-1; d < diffusion.size()-1+Na_i_diffusion.size(); d++) {
								reaction_value = Na_i_diffusion[d-(diffusion.size()-1)];
								reactions[j][i][k][r + d] = reaction_value;

								/*if (reaction_value >= max_reaction) {
									max_reaction = reaction_value;
									cell[0] = j; cell[1] = i; cell[2] = k;

									if (d >= 0 && d <= 2)
										reaction_choice = 8;
									else if (d >= 3 && d <= 5)
										reaction_choice = 9;
									else if (d >= 6 && d <= 8)
										reaction_choice = 10;
									else if (d >= 9 && d <= 11)
										reaction_choice = 11;
									else if (d >= 12 && d <= 14)
										reaction_choice = 12;
									else
										reaction_choice = 13;
								}*/

								alfa_0 += reaction_value;
							}
						}
						

						if ((r != 8) || (r != 9)) {
							/*if (reaction_value >= max_reaction) {
								max_reaction = reaction_value;
								cell[0] = j; cell[1] = i; cell[2] = k;
								reaction_choice = r;
							}*/

							alfa_0 += reaction_value;
						}
					}
					// End Reactions >>
				}
			}
		}


		// for (int i = 0; i < DIM_X; i++) {
		// 	for (int j = 0; j < DIM_Y; j++) {
		// 		for (int k = 0; k < DIM_Z; k++) {
		// 			cout << tecido->getId(i, j, k) << ": ";
		// 			for (int r = 0; r < num_reactions; r++) {
		// 				cout << reactions[i][j][k][r] << " ";
		// 			}
		// 			cout << endl;
		// 		}
		// 	}
		// }


		// Gerando dois números aleatórios: r1 e r2
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator (seed);
		std::uniform_real_distribution<double> distribution (0.0,1.0);

		double r1 = distribution(generator);
		double r2 = distribution(generator);

		// Calculando o tempo tau
		double tau = (1 / alfa_0) * log(1 / r1);

		// Definindo a reação que será executada
		double sum_upper = 0, sum_down = 0;
		bool flag = false;

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				for (int k = 0; k < DIM_Z; k++) {
					for (int r = 0; r < (num_reactions + 2*(nConnections * 3 - 1)); r++) {
						// cout << "R: " << r <<endl;
						sum_upper += reactions[j][i][k][r];

						if (sum_upper >= alfa_0 * r2) {
							//cout << i << " " << j << " " << k << " " << r << endl;

							flag = false;
							sum_down = 0;
							for (int x = 0; x < DIM_X && flag == false; x++) {
								for (int y = 0; y < DIM_Y && flag == false; y++) {
									for (int z = 0; z < DIM_Z && flag == false; z++) {
										for (int n = 0; n < (num_reactions +  2*(nConnections * 3 - 1)) && flag == false; n++) { //0 a 43 -> tem que ir
											// cout << "N: " << r <<endl;

											if (x == i && y == j && z == k && n == r) {
												//cout << x << " " << y << " " << z << " " << n << endl;
												//cout << sum_down << " " << alfa_0 * r2 << " " << sum_upper << endl;
												flag = true;
											} else {
												sum_down += reactions[y][x][z][n];
											}
										}
									}
								}
							}

							if (sum_down < alfa_0 * r2) {
								//cout << sum_upper << " " << alfa_0 * r2 << " " << sum_down << endl;
								retorno[0] = r + 1; // [1, 26]
								retorno[1] = i;
								retorno[2] = j;
								retorno[3] = k;
								retorno[4] = tau;

								return retorno;
							}
						}
					}
				}
			}
		}
	}
};

void simulation(int destination, double frequency, string topologie) {
	Network tecido;
	int tx_x = trunc(DIM_X / 2);
	int tx_y = trunc(DIM_Y / 2);
	int tx_z = trunc(DIM_Z / 2);

	// DEFININDO A TOPOLOGIA DO TECIDO
	if (topologie == "RD") tecido.regularDegree();
	else if (topologie == "LR2") tecido.linkRadius(2);
	else if (topologie == "LR3") tecido.linkRadius(3);

	tecido.set(tx_x, tx_y, tx_z, "C", 0.5*SCALE);
	tecido.set(tx_x, tx_y, tx_z, "Na_i", 20000*SCALE);

	// Print tissue
	tecido.printTissue();

	// Opening file containing the data of calcium concentration that will be plotted
	ofstream exportfile;
	exportfile.open("temp/data.txt");
	tecido.writeFileHeader(exportfile);
	tecido.printTissuef(exportfile, 0); // Writes the tissue's initial state to the file
	

	// Inicializando o Algoritmo de Gillespie
	Gillespie gillespie(&tecido);

	int nConnections = tecido.numberConnections();
	cout << "Connections per cell: " << nConnections << endl;
	vector<double> choice(5);
	vector<int> connections(nConnections), qtd_reactions(8 + 2 * (nConnections * 3));
	double simulation_time = 200, current_time = 0;
	int reaction, int_time = 0;
	bool diffusion_error = false;

	vector<double> C_tx, C_rx;

	while (simulation_time > current_time) {
		choice = gillespie.run();

		// cout << choice[0] << endl;
		// getchar();

		current_time += (choice[4] * 1000);
		reaction = choice[0];

		diffusion_error = false;
		qtd_reactions[reaction - 1]++;


		// DEBUG (To see reaction by reaction)
		/* cout << endl;
		cout << "Reaction: " << reaction << endl;
		cout << "Cell: " << tecido.getId(choice[1], choice[2], choice[3]) << endl;
		tecido.printTissue();
		cout << endl;
		getchar(); */

		//cout << current_time << endl;

		if (trunc(current_time) != int_time) {
			int_time = trunc(current_time);

			cout << "Time: " << int_time << endl;

			// Print Reactions
			cout << "Reactions: ";
			for (int i = 0; i < qtd_reactions.size(); i++) {
				cout << qtd_reactions[i] << " ";
			}
			cout << endl;
			// End Print Reactions

			// Print Tissue
			tecido.printTissue();
			cout << endl;
			// End Print Tissue

			// Write Calcium concentration of tissue on file
			if (DIM_Z == 1) tecido.printTissuef(exportfile, int_time);
			// End Print Tissue

			// Update Gap Junctions
			if (int_time < 50) {
				for (int i = 0; i < DIM_X; i++) {
					for (int j = 0; j < DIM_Y; j++) {
						for (int k = 0; k < DIM_Z; k++) {
							tecido.set(i, j, k, "phl", phl[int_time]);
							tecido.set(i, j, k, "plh", plh[int_time]);
							tecido.set(i, j, k, "phh", phh[int_time]);
						}
					}
				}
			}

			// Calcium oscillations
			if (int_time % ((int) (1 / frequency)) == 0)
				tecido.accumulate(tx_x, tx_y, tx_z, "C", 0.5);
		}

		/* INTRACELULAR REACTIONS */
		if (reaction == 1) {
			tecido.accumulate(choice[1], choice[2], choice[3], "C", ALPHA);
		} else if (reaction == 2) {
			tecido.accumulate(choice[1], choice[2], choice[3], "C", ALPHA);
			tecido.accumulate(choice[1], choice[2], choice[3], "ER", -ALPHA);

			if (tecido.get(choice[1], choice[2], choice[3], "ER") < 0) {
				tecido.changeSignal(choice[1], choice[2], choice[3], "ER");
			}
		} else if (reaction == 3) {
			tecido.accumulate(choice[1], choice[2], choice[3], "C", -ALPHA);
			tecido.accumulate(choice[1], choice[2], choice[3], "ER", ALPHA);

			if (tecido.get(choice[1], choice[2], choice[3], "C") < 0) {
				tecido.changeSignal(choice[1], choice[2], choice[3], "C");
			}
		} else if (reaction == 4) {
			tecido.accumulate(choice[1], choice[2], choice[3], "C", ALPHA);
			tecido.accumulate(choice[1], choice[2], choice[3], "ER", -ALPHA);

			if (tecido.get(choice[1], choice[2], choice[3], "ER") < 0) {
				tecido.changeSignal(choice[1], choice[2], choice[3], "ER");
			}
		} else if (reaction == 5) {
			tecido.accumulate(choice[1], choice[2], choice[3], "C", -ALPHA);
			tecido.accumulate(choice[1], choice[2], choice[3], "ER", ALPHA);

			if (tecido.get(choice[1], choice[2], choice[3], "C") < 0) {
				tecido.changeSignal(choice[1], choice[2], choice[3], "C");
			}
		} else if (reaction == 6) {
			tecido.accumulate(choice[1], choice[2], choice[3], "C", -ALPHA);

			if (tecido.get(choice[1], choice[2], choice[3], "C") < 0) {
				tecido.changeSignal(choice[1], choice[2], choice[3], "C");
			}
		} else if (reaction == 7) {
			tecido.accumulate(choice[1], choice[2], choice[3], "IP3", ALPHA);
		} else if (reaction == 8) {
			tecido.accumulate(choice[1], choice[2], choice[3], "IP3", -ALPHA);

			if (tecido.get(choice[1], choice[2], choice[3], "IP3") < 0) {
				tecido.changeSignal(choice[1], choice[2], choice[3], "IP3");
			}
		}

		/* DIFFUSION REACTIONS */
		// Calcium Diffusion
		else if(reaction >= 9 && reaction < 9 + (nConnections * 3)){
			for (int conn = 0; conn < nConnections; conn++) {
				if (reaction >= 9 + (conn * 3) && reaction <= 11 + (conn * 3)) {
					connections = tecido.getConnections(tecido.getId(choice[1], choice[2], choice[3]));

					if (connections[conn] != -1 && tecido.get(choice[1], choice[2], choice[3], "C") > tecido.get(connections[conn], "C")) {
						tecido.accumulate(choice[1], choice[2], choice[3], "C", -ALPHA);
						tecido.accumulate(connections[conn], "C", ALPHA);
					} else {
						cout << tecido.getId(choice[1], choice[2], choice[3]) << " " << connections[conn] << endl;
						diffusion_error = true;
					}
				}
			}
		}
	
		// Sodium Diffusion
		else {
			for (int conn = 0; conn < nConnections; conn++) {
				if (reaction >= 27 + (conn * 3) && reaction <= 29 + (conn * 3)) {
					connections = tecido.getConnections(tecido.getId(choice[1], choice[2], choice[3]));

					if (connections[conn] != -1 && tecido.get(choice[1], choice[2], choice[3], "Na_i") > tecido.get(connections[conn], "Na_i")) {
						tecido.accumulate(choice[1], choice[2], choice[3], "Na_i", -ALPHA);
						tecido.accumulate(connections[conn], "Na_i", ALPHA);
					} else {
						cout << tecido.getId(choice[1], choice[2], choice[3]) << " " << connections[conn] << endl;
						diffusion_error = true;
					}
				}
			}
		}
		
		/* STORAGE OF CALCIUM CONCENTRATION */
		C_tx.push_back(tecido.get(tx_x, tx_y, tx_z, "C"));
		C_rx.push_back(tecido.get(tx_x + destination, tx_y, tx_z, "C"));

		/* DIFFUSION ERROR => NO INCREMENT TIME */
		if (diffusion_error) current_time -= (choice[4] * 1000);
	}


	/* ### CALCULATING GAIN ### */
	ofstream file_gain;
	file_gain.open("results/gain.txt", ios::app);

	double acc_c_tx = accumulate(C_tx.begin(), C_tx.end(), 0.0);
	double acc_c_rx = accumulate(C_rx.begin(), C_rx.end(), 0.0);

	double calc_gain = 10 * log10((acc_c_rx / C_rx.size()) / ((acc_c_tx) / C_tx.size()));
	file_gain << topologie << " | Rx = " << destination << "; Freq = " << frequency <<  " Hz | Gain = " << calc_gain << endl;

	file_gain.close();

	/* ### END GAIN ### */
}

/* MAIN */
int main(){
	int destination = 1;
	double frequency = 0.000001;
	string topologie = "RD";
	
	//for (int dest = 1; dest <= 1; dest++) {
	//	simulation(dest, frequency, topologie);
	//}

	simulation(destination, frequency, topologie);

	return 0;
}
