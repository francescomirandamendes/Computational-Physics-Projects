#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>      
#include <cstdio>
#include <stdlib.h>
#include <tuple>
#include <math.h>
#include <sstream>

#define ImageWidth 1000  //image width
#define ImageHeight 1000 //image height

using namespace std;

void Print_lattice (int *vlat, const int &vlx, const int &vly, const int &vwidth, const int &vheight, const char* vfilename="output.ppm")
{
  int  i, j, k, l;
  int vw= vwidth/vlx, vh=vheight/vly;
  int r[5], g[5], b[5];

  r[0]= 255; g[0]= 255; b[0]= 255; //white
  r[1]=   0; g[1]= 255; b[1]=   0; //green
  r[2]= 255; g[2]=   0; b[2]=   0; //red
  r[3]=   0; g[3]=   0; b[3]=   0; //black
  r[4]=   0; g[4]=   0; b[4]= 255; //blue

  ofstream out (vfilename);

  out << "P3" << endl;
  out << vw*vlx << " " << vh*vly << endl;
  out << "255" << endl;

  for (i=vly-1; i>=0; i--)
    for (j=0; j<vh; j++)
      for (k=0; k<vlx; k++)
      {
        for (l=0; l<vw; l++)
        { out << r[vlat[k+i*vlx]] << " " << g[vlat[k+i*vlx]] << " " << b[vlat[k+i*vlx]] << " ";
        }
      } 
      out << endl;

  out.close ();
}
//2.1
void config(double p, int seed, int N){ //gerador de configuração
  srand48(seed);
  int lat[N*N];  //create lattice
  for(int i=0; i<N*N; i++){  
     double m = drand48();
      
     if(m<=p){
       lat[i] = 1;}//ocupado
       else{
         lat[i] = 0;}//livre
  }
  Print_lattice (lat, N, N, ImageWidth, ImageHeight, "configuracao.ppm");
}
//2.2
tuple<int, int, bool> queima(double p, int N,  bool save_images = false){
    int lat[N*N];

    // Gerar configuração inicial
    for (int icounter=0; icounter<N*N; icounter++){
        double m = drand48();
        if(m<=p){
            lat[icounter] = 1;  // ocupado
        } else {
            lat[icounter] = 0;  // livre
        }
    }
    //IMPRIMIR FOTOS (NÃO USAR PARA 2.3)
    if(save_images){
    Print_lattice (lat, N, N, ImageWidth, ImageHeight, "burn_00.ppm"); 
    }
    
    // Queima da primeira linha (linha 0)
    for(int icounter=0; icounter<N; icounter++){
        if(lat[icounter]==1){
            lat[icounter]=2;
        }
    }
    
    bool fire = true;
    int time = 0;
    bool cluster = false;
    int path = 0;
    
    while(fire){
        time++;
        
        // Criar array temporário para a próxima iteração
        int lat_new[N*N];
        for(int i=0; i<N*N; i++){
            lat_new[i] = lat[i];
        }
        
        bool burned = false;
        
        // Processar todas as células
        for(int i=0; i<N*N; i++){
            if(lat[i]==2){  // se está a queimar
                burned = true;
                
                // Vizinho da esquerda
                if((i % N) != 0 && lat[i-1]==1){
                    lat_new[i-1]=2;
                }
                
                // Vizinho da direita
                if((i % N) != (N-1) && lat[i+1]==1){
                    lat_new[i+1]=2;
                }
                
                // Vizinho de cima
                if(i+N < N*N && lat[i+N]==1){
                    lat_new[i+N]=2;
                }
                
                // Vizinho de baixo
                if(i-N >= 0 && lat[i-N]==1){
                    lat_new[i-N]=2;
                }
                
                // Marcar como queimado
                lat_new[i]=3;
            }
        }
        
        // Atualizar lattice
        for(int i=0; i<N*N; i++){
            lat[i] = lat_new[i];
        }
        //IMPRIMIR FOTOS (NÃO USAR PARA 2.3)
        if(save_images){
        char filename[50];
        sprintf(filename, "queima_%02d.ppm", time);
        Print_lattice(lat, N, N, ImageWidth, ImageHeight, filename);
        }
        // Verificar se chegou ao topo
        if(!cluster){
            for(int coluna=0; coluna<N; coluna++){
                if(lat[coluna + (N-1)*N] == 2 || lat[coluna + (N-1)*N] == 3){
                    cluster = true;
                    path = time;
                    break;
                }
            }
        }
        
        // Verificar se ainda há fogo
        if(!burned){
            fire = false;
        }
    }
    
    if(!cluster){
        path = 0;
    }
    if(save_images){
        cout << "Guardadas " << time+1 << " imagens (burn_00.ppm até queima_"
             << setfill('0') << setw(2) << time << ".ppm)\n";
        
        if(cluster){
            cout << "Houve percolação, sendo que o percurso mais curto foi de " 
                 << path << " passos.\n";
        } else {
            cout << "Não houve percolação.\n";
        }
    cout << "O fogo parou ao fim de " << time << " iterações.\n";}
        
    return {path, time, cluster};
    

}

tuple<double, double, double> estatisticas(double p, int N, int n_amostras, int seed){
     srand48(seed);
     double short_path = 0;
     double time = 0;
     double num_percola = 0;

     for (int n=0; n<n_amostras; n++){
         auto [x1,x2,x3] = queima(p,N,false);
         time += x2;
         if(x3){
           num_percola++;
           short_path += x1;}
    }
    
     double frac_percola = num_percola / n_amostras;
     double spath_medio  = (num_percola > 0) ? (short_path / num_percola) : 0; 
     double time_medio = time / n_amostras;
     
     cout << "Para N = " << N << ", p = " << p << ", n_amostras = " << n_amostras 
          << " (seed = " << seed << "):\n";
     cout << "  Fração que percola = " << frac_percola << "\n";
     cout << "  Caminho curto médio = " << spath_medio << "\n";
     cout << "  Tempo médio = " << time_medio << "\n\n";
     
  return {frac_percola, time_medio, spath_medio};
     
    }
    
static vector<double> arange(double ini, double fim, double passo) {
    vector<double> v;
    for (double x = ini; x <= fim + 1e-12; x += passo) {
 
        double xr = floor(x * 1000.0 + 0.5) / 1000.0;
        v.push_back(xr);
    }
    return v;
}

void dados_graf(int n_amostras){

const vector<int> x_n = {10, 25, 50, 100};
const vector<double> x_p = arange(0.0, 1.0, 0.01);

int seed_base = 42;
int seed_counter = 0;

for (int N : x_n){
    std::ostringstream n1, n2, n3;
        n1 << "frac_perc_N" << N << ".dat";
        n2 << "tempo_N"     << N << ".dat";
        n3 << "spath_N"     << N << ".dat";
        
        ofstream f1(n1.str()), f2(n2.str()), f3(n3.str());
    for (size_t j=0;j<x_p.size();j++){
    
        int seed = seed_base + seed_counter;
        seed_counter++;
        
        double p = x_p[j];
        auto [frac_perc, tempo_med, spath_med] = estatisticas(p, N, n_amostras, seed);

            f1 << p << ' ' << frac_perc << '\n';
            f2 << p << ' ' << tempo_med << '\n';
            f3 << p << ' ' << spath_med << '\n';
            }
     cout << "Concluído N = " << N << "\n";
     }
}

int main(){
  //config(0.25,1,25); 
  //config(0.5,1,25); 
 //config(0.75,1,25);
  queima(0.5,50,true);
  queima(0.7,50,true);
  //dados_graf(1000);
  //estatisticas(0.6, 100,1000,33);
  
  return 0;}
