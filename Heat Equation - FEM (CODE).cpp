#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "sparselib.hh"

using namespace std;
using namespace ::sparselib_load;
using ::sparselib::index_type;


double solucao_analitica(double x, double k, double L, double T_A) {
    double sqrt_k = sqrt(k);
    double exp_2kL = exp(2.0 * sqrt_k * L);
    return T_A * (exp(sqrt_k * (2*L - x)) - exp(sqrt_k * x)) / (exp_2kL - 1.0);
}

// Monta a matriz A e vetor b
void montar_sistema(CCoorMatrix<double>& A, Vector<double>& b, 
                    int N, double delta, double k, double T_A) {
    b = 0.0;
    b[0] = T_A * (1.0/delta - k*delta/6.0);

    for (int i = 0; i < N-1; i++)
        A.insert(i, i) = 2.0/delta + 2.0*k*delta/3.0;

    A.insert(0, 1) = -1.0/delta + k*delta/6.0;
    A.insert(N-2, N-3) = -1.0/delta + k*delta/6.0;

    for (int i = 1; i < N-2; i++) {
        A.insert(i, i-1) = -1.0/delta + k*delta/6.0;
        A.insert(i, i+1) = -1.0/delta + k*delta/6.0;
    }
    
    A.internal_order();
}

void calcular_erros(Vector<double>& x, int N, double delta, 
                   double k, double L, double T_A,
                   double& max_erro_abs, double& max_erro_rel) {
    max_erro_abs = 0.0;
    max_erro_rel = 0.0;
    
    for (int i = 0; i < N-1; i++) {
        double x_pos = (i+1) * delta;
        double u_ana = solucao_analitica(x_pos, k, L, T_A);
        double erro_abs = fabs(x[i] - u_ana);
        double erro_rel = fabs(erro_abs / u_ana);
        
        max_erro_abs = max(max_erro_abs, erro_abs);
        max_erro_rel = max(max_erro_rel, erro_rel);
    }
}

void salvar_dados(const string& filename, Vector<double>& x, 
                  int N, double delta, double k, double L, double T_A) {
    ofstream arquivo(filename);
    arquivo << "# x u_numerica u_analitica erro_absoluto erro_relativo" << endl;
    arquivo << setprecision(10);
    
    // Ponto inicial (fronteira)
    double u_ana_0 = solucao_analitica(0.0, k, L, T_A);
    arquivo << 0.0 << " " << T_A << " " << u_ana_0 << " " 
            << fabs(T_A - u_ana_0) << " " << fabs((T_A - u_ana_0)/u_ana_0) << endl;

    for (int i = 0; i < N-1; i++) {
        double x_pos = (i+1) * delta;
        double u_ana = solucao_analitica(x_pos, k, L, T_A);
        double erro_abs = fabs(x[i] - u_ana);
        double erro_rel = fabs(erro_abs / u_ana);
        
        arquivo << x_pos << " " << x[i] << " " << u_ana << " " 
                << erro_abs << " " << erro_rel << endl;
    }
    
    arquivo.close();
}

void resolver_sistema(int N, double epsilon, double k, double L, double T_A,
                     int max_iter, double& max_erro_abs, double& max_erro_rel,
                     int& iter, double& res, Vector<double>& x) {
    double delta = L / N;
    
    Vector<double> b;
    CCoorMatrix<double> A;
    
    x.new_dim(N-1);
    b.new_dim(N-1);
    A.new_dim(N-1, N-1, N-1 + 2*(N-2));
    
    montar_sistema(A, b, N, delta, k, T_A);
    
    x = 0.0;
    DPreco<double> P;
    P.build(A);
    
    res = bicgstab(A, b, x, P, epsilon, max_iter, iter);
    
    calcular_erros(x, N, delta, k, L, T_A, max_erro_abs, max_erro_rel);
}

int main() {
   
    double L = 10.0;
    double k = 1.0;
    double T_A = 1.0;
    double T_R = 2.0 * T_A;
    int max_iter = 1e4;
    double epsi_padrao = 1e-8;

    cout << "Parâmetros:" << endl;
    cout << "  L = " << L << endl;
    cout << "  k = " << k << endl;
    cout << "  T_A = " << T_A << endl;
    cout << "  T_R = " << T_R << endl;
    

    cout << "PARTE 1: Análise para diferentes N (ε = " << epsi_padrao << ")\n" << endl;
    
    int N_values[] = {10, 25, 50, 100};
    int num_N = 4;
    
    ofstream stats("estatisticas.dat");
    stats << "# N erro_max_abs erro_max_rel iteracoes residuo" << endl;
    
    for (int idx = 0; idx < num_N; idx++) {
        int N = N_values[idx];
        double delta = L / N;
        
        cout << "N = " << N << " (Δx = " << delta << "): ";
        
        Vector<double> x;
        double max_erro_abs, max_erro_rel;
        int iter;
        double res;
        
        resolver_sistema(N, epsi_padrao, k, L, T_A, max_iter, 
                        max_erro_abs, max_erro_rel, iter, res, x);
        
        cout << iter << " iter, res = " << ", erro = " << max_erro_abs << fixed << endl;
        
        stats << N << " " << max_erro_abs << " " << max_erro_rel << " " 
              << iter << " "  << endl;
        
        salvar_dados("dados_N" + to_string(N) + ".dat", x, N, delta, k, L, T_A);
    }
    
    stats.close();
    cout << "\n✓ Dados de diferentes N salvos\n" << endl;

    double k_values[] = {0.1, 0.5, 1.0, 2.0, 5.0};
int num_k = 5;
int N_k = 100;
double delta_k = L / N_k;
ofstream dados_k("dados_diferentes_k.dat");
dados_k << "# x k u_numerica u_analitica" << endl;
dados_k << setprecision(10);

ofstream stats_k("estatisticas_k.dat");
stats_k << "# k erro_max_abs erro_max_rel iteracoes residuo" << endl;

for (int idx = 0; idx < num_k; idx++) {
    double k_atual = k_values[idx];
    
    cout << "k = " << k_atual << ": ";
    
    Vector<double> x;
    double max_erro_abs, max_erro_rel;
    int iter;
    double res;
    resolver_sistema(N_k, epsi_padrao, k_atual, L, T_A, max_iter,
                    max_erro_abs, max_erro_rel, iter, res, x);
    
    cout << iter << " iter, erro = " << scientific << max_erro_abs << fixed << endl;
    
    stats_k << k_atual << " " << max_erro_abs << " " << max_erro_rel << " " 
            << iter << " " << res << endl;
    
    dados_k << 0.0 << " " << k_atual << " " << T_A << " " 
            << solucao_analitica(0.0, k_atual, L, T_A) << endl;

    for (int i = 0; i < N_k-1; i++) {
        double x_pos = (i+1) * delta_k;
        dados_k << x_pos << " " << k_atual << " " << x[i] << " " 
                << solucao_analitica(x_pos, k_atual, L, T_A) << endl;
    }

    dados_k << L << " " << k_atual << " " << 0.0 << " " 
            << solucao_analitica(L, k_atual, L, T_A) << endl;

    salvar_dados("dados_k" + to_string(k_atual) + ".dat", x, N_k, delta_k, k_atual, L, T_A);
}

dados_k.close();
stats_k.close();
cout << "\n✓ Dados de diferentes k salvos" << endl;

    
    cout << "  Arquivos gerados:" << endl;
    cout << "  - dados_N*.dat (4 arquivos)" << endl;
    cout << "  - dados_diferentes_epsilon.dat" << endl;
    cout << "  - estatisticas.dat" << endl;
    cout << "  - estatisticas_epsilon.dat" << endl;
    cout << "========================================" << endl;
    
    return 0;
}

