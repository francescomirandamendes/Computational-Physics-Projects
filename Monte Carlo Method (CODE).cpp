#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <string>
#include <bits/stdc++.h>
#include <sstream>
using namespace std;

double generator(long double mu_tot){
      long double x_i = -(1/mu_tot) * log(1-drand48());
      return (x_i);
      }
      
void fotoes1(){
    
    cout << "EXERCICIO 1: Monocromatico" <<" \n";
     
    double mu_f = 3.19 * pow(10,3);
    double mu_R = 6.09;
    double mu_C = 0.0383;
    double mu_tot = mu_C + mu_f + mu_R;
     
    cout << "Coeficiente total: mu_tot = " << mu_tot << " cm^-1\n";
    cout << "Comprimento de atenuação: lambda = " << 1.0/mu_tot << " cm\n\n";
    
    double I_0 = 100000;
    double n_bins = 200;
    vector<double> x_memory;
    
    for (long i = 0; i<I_0; i++){
        long double xi = generator(mu_tot);
        x_memory.push_back(xi);
        }
        
        
    sort(x_memory.begin(), x_memory.end());
    double x_max = x_memory[x_memory.size()-1];
    cout << "O valor máximo de x é "<< x_max <<"\n";
    double delta_x = x_max/n_bins;
    cout << "A largura dos bins é de "<<delta_x<<" cm"<<"\n";
    
    
    ofstream file1("dados_histograma_1.txt");
    ofstream file2("dados_histograma_1log.txt");
    for (int bin = 0; bin < n_bins; bin++){
        double x = (bin + 0.5) * delta_x;
        int sobreviventes = 0;
        
        for (int j = 0; j<I_0; j++){
            if (x_memory[j] > x){
               sobreviventes++;}}
         
        double hist_x = (double)sobreviventes / I_0;
        double logx = log((double)sobreviventes / I_0);
        file1 << x << "\t"<< hist_x << "\n";
        file2 << x << "\t" <<logx<<"\n";
    }
    
    file1.close();
    file2.close();
    vector<double> x_sorted = x_memory;
    sort(x_sorted.begin(), x_sorted.end());
    double hvl = x_sorted[I_0/2 - 1];
   cout<<"O ponto x para o qual I/I_0 é metade é: "<<hvl<<"\n";
       }

void fotoes2(){
    ifstream fileSpec("110kV.spec");
    cout << "EXERCICIO 2: Policromatico" <<" \n";
    vector<double> energiasSpec;
    vector<double> fluxoSpec;
    string linha;
    bool lendoSpec = false;
    
    while(getline(fileSpec, linha)){
        if(linha.find("**** CALCULATED SPECTRUM ****") != string::npos){
          lendoSpec = true;
          getline(fileSpec, linha);
          continue;
         }
         
     if(lendoSpec && !linha.empty() && linha[0] != '*'){
       istringstream iss(linha);
       double E,N;
       if(iss >> E >> N){
         energiasSpec.push_back(E);
         fluxoSpec.push_back(N);
         }
       }
    }
    fileSpec.close();
    
       
       
    ifstream fileGa("Galium.txt");
    vector<double> energiasGa; //MeV
    vector<double> muR_rho;
    vector<double> muC_rho;
    vector<double> muF_rho;
    
    double densidade = 5.904; // g/cm^3
    
    getline(fileGa, linha);
    getline(fileGa, linha);
    
    while(getline(fileGa, linha)){
        if(linha.empty()) continue;
        
        istringstream iss(linha);
        double E, muR, muC, muF;
        
        if (iss>>E>>muR>>muC>>muF){
           energiasGa.push_back(E*1000);
           muR_rho.push_back(muR);
           muC_rho.push_back(muC);
           muF_rho.push_back(muF);
           }
     }
    fileGa.close();
    
    vector<double> mu_R(energiasGa.size());
    vector<double> mu_C(energiasGa.size());
    vector<double> mu_F(energiasGa.size());
    vector<double> mu_tot(energiasGa.size());
    
    for(int i = 0; i<energiasGa.size(); i++){
       mu_R[i] = muR_rho[i] * densidade;
       mu_C[i] = muC_rho[i] * densidade;
       mu_F[i] = muF_rho[i] * densidade;
       mu_tot[i] = mu_R[i] + mu_C[i] + mu_F[i];
       }
     vector<double> probR(energiasGa.size());
     vector<double> probC(energiasGa.size());
     vector<double> probF(energiasGa.size());
     
     for (int i = 0; i < energiasGa.size(); i++){
        probR[i] = mu_R[i]/mu_tot[i];
        probC[i] = mu_C[i]/mu_tot[i];
        probF[i] = mu_F[i]/mu_tot[i];
        }
     double I_0 = 10000;

     vector<double> x_memory;
     vector<int> processo_memory; //0=R 1=C 2=F
     
     double somaFluxo = 0;
    for(size_t i = 0; i < fluxoSpec.size(); i++){
        somaFluxo += fluxoSpec[i];
    }

    vector<int> nFotoes_por_energia(energiasSpec.size());
    
    for(size_t i = 0; i < energiasSpec.size(); i++){
        // Proporção de fotões com esta energia
        double frac = fluxoSpec[i] / somaFluxo;
        nFotoes_por_energia[i] = (int)(frac * I_0);
    }
    
    int soma_parcial = 0;
    for(size_t i = 0; i < nFotoes_por_energia.size(); i++){
        soma_parcial += nFotoes_por_energia[i];
    }
    nFotoes_por_energia[0] += (I_0 - soma_parcial); // Ajuste no primeiro bin
    
    int contador = 0;

    for(size_t iE = 0; iE < energiasSpec.size(); iE++){
        double E = energiasSpec[iE];
        int nFotoes = nFotoes_por_energia[iE];
        
        if(nFotoes == 0) continue;
        double mu_tot_E, pR_E, pC_E, pF_E;
        
        // Procura por energia mais próxima
        double diff_min = 1e10;
        size_t idx_min = 0;
        for(size_t j = 0; j < energiasGa.size(); j++){
            double diff = abs(energiasGa[j] - E);
            if(diff < diff_min){
                diff_min = diff;
                idx_min = j;
            }
        }
        
        mu_tot_E = mu_tot[idx_min];
        pR_E = probR[idx_min];
        pC_E = probC[idx_min];
        pF_E = probF[idx_min];
        
        // Gera os fotões com esta energia
        for(int i = 0; i < nFotoes; i++){
            double x = generator(mu_tot_E);
            x_memory.push_back(x);
            
            // Determina processo
            double r = drand48();
            if(r < pF_E){
                processo_memory.push_back(2); // Fotoeléctrico
            } else if(r < pF_E + pC_E){
                processo_memory.push_back(1); // Compton
            } else {
                processo_memory.push_back(0); // Rayleigh
            }
            
            contador++;
        }
    }
    
    vector<double> x_sorted = x_memory;
    sort(x_sorted.begin(), x_sorted.end());
    double hvl = x_sorted[contador/2 - 1];
    
    cout << "   HVL (espessura para 50% de transmissão): " << hvl << " cm\n";
    
    int n_bins = 200;
    double x_max = *max_element(x_memory.begin(), x_memory.end());
    double delta_x = x_max / n_bins;
    
    cout << "   Máximo de x: " << x_max << " cm\n";
    cout << "   Largura dos bins: " << delta_x << " cm\n";

    vector<int> bins(n_bins, 0);
    vector<int> binsR(n_bins, 0);
    vector<int> binsC(n_bins, 0);
    vector<int> binsF(n_bins, 0);
    
    for(size_t i = 0; i < x_memory.size(); i++){
        int bin = int(x_memory[i] / delta_x);
        if(bin >= n_bins) bin = n_bins - 1;
        
        bins[bin]++;
        
        if(processo_memory[i] == 0) binsR[bin]++;
        if(processo_memory[i] == 1) binsC[bin]++;
        if(processo_memory[i] == 2) binsF[bin]++;
    }
    
    // Acumulação
    vector<int> sobram(n_bins);
    vector<int> acumR(n_bins);
    vector<int> acumC(n_bins);
    vector<int> acumF(n_bins);
    
    sobram[0] = contador - bins[0];
    acumR[0] = binsR[0];
    acumC[0] = binsC[0];
    acumF[0] = binsF[0];
    
    for(int i = 1; i < n_bins; i++){
        sobram[i] = sobram[i-1] - bins[i];
        acumR[i] = acumR[i-1] + binsR[i];
        acumC[i] = acumC[i-1] + binsC[i];
        acumF[i] = acumF[i-1] + binsF[i];
    }

    ofstream file1("dados_histograma_2.txt");
    ofstream file2("dados_processos_2.txt");
    ofstream file3("dados_histograma_2log.dat");
    
    for(int i = 0; i < n_bins - 1; i++){
        double x = (i + 0.5) * delta_x;
        double fracao = (double)sobram[i] / contador;
        
        file1 << x << "\t" << fracao << "\n";
        
        if(fracao > 0){
            file3 << x << "\t" << log(fracao) << "\n";
        }
        
        file2 << x << "\t" 
              << (double)acumR[i]/contador << "\t"
              << (double)acumC[i]/contador << "\t"
              << (double)acumF[i]/contador << "\n";
    }
    
    file1.close();
    file2.close();
    file3.close();
    
    double soma_x = 0;
    for(size_t i = 0; i < x_memory.size(); i++){
        soma_x += x_memory[i];
    }
    double lambda_medio = soma_x / contador;
    cout << "Lambda médio: " << lambda_medio << " cm\n";
    cout << "Mu médio:     " << 1.0/lambda_medio << " cm^-1\n";
    }   
     
int main(){
   srand48(6);
   fotoes1();
   fotoes2();
   
   return 0;
   }
       
       
