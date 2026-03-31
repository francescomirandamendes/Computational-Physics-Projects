#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <set>
  
using namespace std;

// GERADOR 1: Linear Congruential Generator
void gerar_LCG(long int c, int a, long int p, int x0, int N, const string& file2D, const string& file3D) {
    ofstream f2(file2D), f3(file3D);
    f2 << fixed << setprecision(7);
    f3 << fixed << setprecision(7);
    long int x = x0;
    
    set<long int> valores_unicos;
    
    for (int i = 0; i < N; ++i) {
        long long int x1 = (c * x + a) % p;
        long long int x2 = (c * x1 + a) % p;
        double u = static_cast<double>(x) / p;
        double v = static_cast<double>(x1) / p;
        double z = static_cast<double>(x2) / p;
        
        f2 << u << ' ' << v << '\n';
        f3 << u << ' ' << v << ' ' << z << '\n';
        
        valores_unicos.insert(x);
        x = x1;
    }
    
    cout << "  Números únicos gerados: " << valores_unicos.size() << " de " << N << endl;
    cout << "  Arquivos " << file2D << " e " << file3D << " criados." << endl;
}

// GERADOR 2: rand()
void gerar_rand(int N, int seed, const string& file2D, const string& file3D) {
    ofstream f2(file2D), f3(file3D);
    srand(seed);
    f2 << fixed << setprecision(7);
    f3 << fixed << setprecision(7);
    
    set<int> valores_unicos;
    
    int x0_int = rand();
    double x0 = static_cast<double>(x0_int) / RAND_MAX;
    valores_unicos.insert(x0_int);
    
    for (int i = 0; i < N; ++i) {
        int x1_int = rand();
        int x2_int = rand();
        double x1 = static_cast<double>(x1_int) / RAND_MAX;
        double x2 = static_cast<double>(x2_int) / RAND_MAX;
        f2 << x0 << ' ' << x1 << '\n';
        f3 << x0 << ' ' << x1 << ' ' << x2 << '\n';
        valores_unicos.insert(x1_int);
        x0 = x1;
    }
    
    cout << "  RAND_MAX = " << RAND_MAX << endl;
    cout << "  Números únicos gerados: " << valores_unicos.size() << " de " << (N+1) << endl;
    cout << "  Arquivos " << file2D << " e " << file3D << " criados." << endl;
}

// GERADOR 3: drand48()
void gerar_drand48(int N, int seed, const string& file2D, const string& file3D) {
    ofstream f2(file2D), f3(file3D);
    srand48(seed);
    f2 << fixed << setprecision(7);
    f3 << fixed << setprecision(7);
    
    // drand48 usa 48 bits internamente, vamos contar valores double únicos
    set<double> valores_unicos;
    
    double x0 = drand48();
    valores_unicos.insert(x0);
    
    for (int i = 0; i < N; ++i) {
        double x1 = drand48();
        double x2 = drand48();
        f2 << x0 << ' ' << x1 << '\n';
        f3 << x0 << ' ' << x1 << ' ' << x2 << '\n';
        valores_unicos.insert(x1);
        x0 = x1;
    }
    
    cout << "  Números únicos gerados: " << valores_unicos.size() << " de " << (N+1) << endl;
    cout << "  Arquivos " << file2D << " e " << file3D << " criados." << endl;
}

// PONTOS NUM CÍRCULO
void circulo_pts(int c, int a, int p, int x0r, int x0theta, int N, const string& fileR, const string& fileSQRT) {
    ofstream f1(fileR), f2(fileSQRT);
    f1 << fixed << setprecision(7);
    f2 << fixed << setprecision(7);
    
    for (int i = 0; i < N; ++i) {
        long long int x1 = (c * x0r + a) % p;
        
        double xR = static_cast<double>(x1) / p;
        double xSQRT = sqrt(static_cast<double>(x1) / p);
        
        long long int x2 = (c * x0theta + a) % p;
        double theta = 2 * M_PI * (static_cast<double>(x2) / p);
        
        double Xr = xR * cos(theta);
        double Yr = xR * sin(theta);
        
        double Xsqrt = xSQRT * cos(theta);
        double Ysqrt = xSQRT * sin(theta);
        
        f1 << Xr << ' ' << Yr << '\n';
        f2 << Xsqrt << ' ' << Ysqrt << '\n';
        
        x0r = x1;
        x0theta = x2;
    }
    cout << "  Arquivos " << fileR << " e " << fileSQRT << " criados." << endl;
}

// TESTE DE CHI-QUADRADO (Exercício 1.3)
void teste_chi2_LCG(long int c, int a, long int p, int x0, int N, int k, const string& outputFile) {
    vector<int> freq(k, 0);
    double esperado = static_cast<double>(N) / k;
    
    long int x = x0;
    
    // Gerar números e contar frequências
    for (int i = 0; i < N; ++i) {
        x = (c * x + a) % p;
        double u = static_cast<double>(x) / p;
        
        int intervalo = static_cast<int>(u * k);
        if (intervalo == k) intervalo = k - 1; // Caso u == 1.0
        freq[intervalo]++;
    }
    
    // Calcular chi-quadrado
    double chi2 = 0.0;
    for (int i = 0; i < k; ++i) {
        double diff = freq[i] - esperado;
        chi2 += (diff * diff) / esperado;
    }
    
    // Mandar para .DAT
    ofstream fout(outputFile);
    fout << "# Teste de Chi-Quadrado\n";
    fout << "# N = " << N << ", k = " << k << "\n";
    fout << "# Chi² = " << fixed << setprecision(6) << chi2 << "\n";
    fout << "# Graus de liberdade = " << (k-1) << "\n";
    fout << "# Frequência esperada = " << esperado << "\n";
    fout << "#\n";
    fout << "# Intervalo    Observado    Esperado    Contribuição\n";
    
    for (int i = 0; i < k; ++i) {
        double contrib = pow(freq[i] - esperado, 2) / esperado;
        fout << setw(6) << i << "    "
             << setw(10) << freq[i] << "    "
             << setw(10) << fixed << setprecision(2) << esperado << "    "
             << setw(12) << setprecision(6) << contrib << "\n";
    }
    fout.close();
    
    // Imprimir terminal
    
    cout << "  Chi² calculado = " << fixed << setprecision(6) << chi2 << endl;
    cout << "  Graus de liberdade = " << (k-1) << endl;
    cout << "  Resultados guardados em: " << outputFile << endl;
}

// Teste chi² para rand()
void teste_chi2_rand(int N, int seed, int k, const string& outputFile) {
    vector<int> freq(k, 0);
    double esperado = static_cast<double>(N) / k;
    
    srand(seed);
    
    for (int i = 0; i < N; ++i) {
        double u = static_cast<double>(rand()) / RAND_MAX;
        int intervalo = static_cast<int>(u * k);
        if (intervalo == k) intervalo = k - 1;
        freq[intervalo]++;
    }
    
    double chi2 = 0.0;
    for (int i = 0; i < k; ++i) {
        double diff = freq[i] - esperado;
        chi2 += (diff * diff) / esperado;
    }
    
    ofstream fout(outputFile);
    fout << " Teste de Chi-Quadrado (rand())\n";
    fout << " N = " << N << ", k = " << k << "\n";
    fout << " Chi² = " << fixed << setprecision(6) << chi2 << "\n";
    fout << " Graus de liberdade = " << (k-1) << "\n";
    fout << "#\n";
    fout << " Intervalo    Observado    Esperado    Contribuição\n";
    
    for (int i = 0; i < k; ++i) {
        double contrib = pow(freq[i] - esperado, 2) / esperado;
        fout << setw(6) << i << "    "
             << setw(10) << freq[i] << "    "
             << setw(10) << fixed << setprecision(2) << esperado << "    "
             << setw(12) << setprecision(6) << contrib << "\n";
    }
    fout.close();
    
  
    cout << "  Chi² calculado = " << fixed << setprecision(6) << chi2 << endl;
    cout << "  Resultados guardados em: " << outputFile << endl;
}

// Teste chi² para drand48()
void teste_chi2_drand48(int N, int seed, int k, const string& outputFile) {
    vector<int> freq(k, 0);
    double esperado = static_cast<double>(N) / k;
    
    srand48(seed);
    
    for (int i = 0; i < N; ++i) {
        double u = drand48();
        int intervalo = static_cast<int>(u * k);
        if (intervalo == k) intervalo = k - 1;
        freq[intervalo]++;
    }
    
    double chi2 = 0.0;
    for (int i = 0; i < k; ++i) {
        double diff = freq[i] - esperado;
        chi2 += (diff * diff) / esperado;
    }
    
    ofstream fout(outputFile);
    fout << " Teste de Chi-Quadrado (drand48())\n";
    fout << " N = " << N << ", k = " << k << "\n";
    fout << " Chi² = " << fixed << setprecision(6) << chi2 << "\n";
    fout << " Graus de liberdade = " << (k-1) << "\n";
    fout << "\n";
    fout << " Intervalo    Observado    Esperado    Contribuição\n";
    
    for (int i = 0; i < k; ++i) {
        double contrib = pow(freq[i] - esperado, 2) / esperado;
        fout << setw(6) << i << "    "
             << setw(10) << freq[i] << "    "
             << setw(10) << fixed << setprecision(2) << esperado << "    "
             << setw(12) << setprecision(6) << contrib << "\n";
    }
    fout.close();
    
    
    cout << "  Chi² calculado = " << fixed << setprecision(6) << chi2 << endl;
    cout << "  Resultados guardados em: " << outputFile << endl;
}

// FUNÇÃO PRINCIPAL
int main() {
    cout << "\n" << endl;
    cout << "  EXERCÍCIO 1: GERADORES ALEATÓRIOS" << endl;
    cout << "\n" << endl;
    
    // Parâmetros
    int c = 3;
    int c_bom = 16807;
    int a = 2;
    int p = 31;
    int p_bom = 2147483647;
    int x0 = 5;
    int seed = 5;
    
    // Exercício 1.1: Teste do Quadrado e do Cubo
    cout << "EXERCÍCIO 1.1: Testes de Correlação\n" << endl;
    
    cout << "1. LCG com MAUS parâmetros (c=" << c << ", p=" << p << "):" << endl;
    int N1 = 1000;
    gerar_LCG(c, a, p, x0, N1, "LCG_2D_mau.dat", "LCG_3D_mau.dat");
    
    cout << "\n2. LCG com BONS parâmetros (c=" << c_bom << ", p=" << p_bom << "):" << endl;
    gerar_LCG(c_bom, a, p_bom, x0, N1, "LCG_2D_bom.dat", "LCG_3D_bom.dat");
    
    cout << "\n3. Gerador rand() C++:" << endl;
    int N2 = 1000;
    gerar_rand(N2, seed, "rand_2D.dat", "rand_3D.dat");
    
    cout << "\n4. Gerador drand48():" << endl;
    gerar_drand48(N2, seed, "drand48_2D.dat", "drand48_3D.dat");
    
    // Exercício 1.2: Círculo
    cout << "\n\nEXERCÍCIO 1.2: Pontos Uniformes no Círculo\n" << endl;
    int N3 = 5000;
    int x0r = 145;
    int x0theta = 698;
    circulo_pts(c_bom, a, p_bom, x0r, x0theta, N3, "circulo_r.dat", "circulo_sqrt.dat");
    
    // Exercício 1.3: Teste de Chi-Quadrado
    cout << "\n\nEXERCÍCIO 1.3: Teste de χ²\n" << endl;
    
    int N_chi = 10000;
    int k = 21;
    int k2 = 51;
    
    cout << " TESTE χ² PARA LCG (MAUS PARÂMETROS):" << endl;
    cout << "K = 21: " << endl;
    teste_chi2_LCG(c, a, p, x0, N_chi, k, "chi2_LCG_mau_k21.dat");
    cout << "K = 51: " << endl;
    teste_chi2_LCG(c, a, p, x0, N_chi, k2, "chi2_LCG_mau_k51.dat");
    
    
    cout << "\n TESTE χ² PARA LCG (BONS PARÂMETROS):" << endl;
    cout << "K = 21: " << endl;
    teste_chi2_LCG(c_bom, a, p_bom, x0, N_chi, k, "chi2_LCG_bom_k21.dat");
    cout << "K = 51: " << endl;
    teste_chi2_LCG(c_bom, a, p_bom, x0, N_chi, k2, "chi2_LCG_bom_k51.dat");
    
    cout << "\n TESTE χ² PARA RAND():" << endl;
    cout << "K = 21: " << endl;
    teste_chi2_rand(N_chi, seed, k, "chi2_rand_k21.dat");
    cout << "K = 51: " << endl;
    teste_chi2_rand(N_chi, seed, k2, "chi2_rand_k51.dat");
    
    cout << "\n TESTE χ² PARA DRAND48():" << endl;
    cout << "K = 21: " << endl;
    teste_chi2_drand48(N_chi, seed, k, "chi2_drand48_k21.dat");
    cout << "K = 51: " << endl;
    teste_chi2_drand48(N_chi, seed, k2, "chi2_drand48_k51.dat");
    
    cout << "\n" << endl;
    cout << "  TODOS OS FICHEIROS FORAM CRIADOS!" << endl;
    cout << "\n" << endl;
    
    return 0;
}
