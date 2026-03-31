#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>
#include <sstream>

using namespace std;

const double G = 4.4939e-15;
const double VELCONV = 1.0226911586616e-6;

inline double distance_sq(double dx, double dy) {
    return dx*dx + dy*dy;
}

inline double distance(double dx, double dy) {
    return sqrt(dx*dx + dy*dy);
}

double gerar_massa(double m_min, double m_max, double alpha){
    double u = drand48();
    double exp_alpha = 1.0 - alpha;
    double m_min_pow = pow(m_min, exp_alpha);
    double m_max_pow = pow(m_max, exp_alpha);
    double m = pow(m_min_pow + u * (m_max_pow - m_min_pow), 1.0/exp_alpha);
    return m;
}
    
double perfil_king(double r, double n0, double rc, double rt) {
    double term1 = pow(1.0 + (r/rc)*(r/rc), -0.5);
    double term2 = pow(1.0 + (rt/rc)*(rt/rc), -0.5);
    return n0 * (term1 - term2);
}

double gerar_raio_king(double rc, double rt) {
    double n0 = 1.0;
    double n_max = perfil_king(0.0, n0, rc, rt);
    
    while (true) {
        double r = drand48() * rt;
        double u = drand48();
        double n_r = perfil_king(r, n0, rc, rt);
        
        if (u * n_max < n_r) {
            return r;
        }
    }
}

double randn() {
    double u1 = drand48();
    double u2 = drand48();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}
    
void um_ponto_um() {
    const int N = 2;
    const int iter = 1800;
    const double dt = 1e5;
   
    vector<double> x_old(N);
    vector<double> y_old(N);
    vector<double> x_now(N);
    vector<double> y_now(N);
    
    x_old[0] = -1.0;  x_old[1] = 1.0;  
    y_old[0] = 0.0;   y_old[1] = 0.0;
    
    double vx0 = 0.0, vx1 = 0.0;
    double vy0 = 0.02 * VELCONV;
    double vy1 = -0.02 * VELCONV;
    
    double m0 = 1.0, m1 = 1.0;
    
    x_now[0] = x_old[0] + vx0 * dt;
    x_now[1] = x_old[1] + vx1 * dt;
    y_now[0] = y_old[0] + vy0 * dt;
    y_now[1] = y_old[1] + vy1 * dt;
    
    vector<double> E_k(iter);
    vector<double> E_p(iter);
    vector<double> E_m(iter);
    vector<double> t_milhoes(iter);

    ofstream positions("ex1_pos.txt");
    ofstream energies("ex1_ngr.txt");
    
    for (int i = 0; i < iter; i++) {
        double t_Myr = i * dt * 1e-6;
        t_milhoes[i] = t_Myr;
  
        double dx = x_now[1] - x_now[0];
        double dy = y_now[1] - y_now[0];
        double r = distance(dx, dy);
        
        double r_inv = 1.0 / r;
        double r_inv3 = r_inv * r_inv * r_inv;
        
        double ax_factor = G * dx * r_inv3;
        double ay_factor = G * dy * r_inv3;
        
        double ax0 = m1 * ax_factor;
        double ay0 = m1 * ay_factor;
        double ax1 = -m0 * ax_factor;
        double ay1 = -m0 * ay_factor;
        
        double x_next0 = 2.0 * x_now[0] - x_old[0] + ax0 * dt * dt;
        double y_next0 = 2.0 * y_now[0] - y_old[0] + ay0 * dt * dt;
        double x_next1 = 2.0 * x_now[1] - x_old[1] + ax1 * dt * dt;
        double y_next1 = 2.0 * y_now[1] - y_old[1] + ay1 * dt * dt;
        
        double dt_inv = 1.0 / (2.0 * dt);
        double vx_0 = (x_next0 - x_old[0]) * dt_inv;
        double vy_0 = (y_next0 - y_old[0]) * dt_inv;
        double vx_1 = (x_next1 - x_old[1]) * dt_inv;
        double vy_1 = (y_next1 - y_old[1]) * dt_inv;
        
        double v0_sq = vx_0*vx_0 + vy_0*vy_0;
        double v1_sq = vx_1*vx_1 + vy_1*vy_1;
        
        double Ec = 0.5 * (m0 * v0_sq + m1 * v1_sq);
        double Ep = -G * m0 * m1 * r_inv;
        double Em = Ec + Ep;
        
        E_k[i] = Ec;
        E_p[i] = Ep;
        E_m[i] = Em;
        
        positions << x_next0 << "\t" << y_next0 << "\t" 
                  << x_next1 << "\t" << y_next1 << endl;
        
        x_old[0] = x_now[0];  x_old[1] = x_now[1];
        y_old[0] = y_now[0];  y_old[1] = y_now[1];
        x_now[0] = x_next0;   x_now[1] = x_next1;
        y_now[0] = y_next0;   y_now[1] = y_next1;
    }
    
    for (int i = 0; i < iter; i++) {
        energies << t_milhoes[i] << "\t" << E_k[i] << "\t" << E_p[i] << "\t" << E_m[i] << endl;
    }
    
    positions.close();
    energies.close();
    
    double E_inicial = E_m[0];
    double E_final = E_m[iter-1];
    double variacao = fabs((E_final - E_inicial) / E_inicial) * 100;
    
    cout << "Energia inicial: " << E_inicial << endl;
    cout << "Energia final: " << E_final << endl;
    cout << "Variação: " << variacao << "%" << endl;
}

void ex2_from_file(const string& input_file, const string& output_file) {
    ifstream data(input_file);
    
    vector<double> x, y, z;           
    vector<double> x_old, y_old, z_old;
    vector<double> ax, ay, az;        
    vector<double> massa;            
    
    string line;
    double total_mass = 0.0;
    
    const double dt = 1e4;
    const double dt_sq = dt * dt;
    const int    iter = 1000;
    const double r_t = 5.0;
    const double r_t_sq = r_t * r_t;
    
    while (getline(data, line)) {
        if (line.size() == 0 || line[0] == '#') continue;
        istringstream iss(line);
        double x_temp, y_temp, z_temp;
        double vx, vy, vz, m;
        
        if (!(iss >> x_temp >> y_temp >> z_temp >> vx >> vy >> vz >> m)) continue;
        
        x_old.push_back(x_temp);
        y_old.push_back(y_temp);
        z_old.push_back(z_temp);
        
        x.push_back(x_temp + vx * dt);
        y.push_back(y_temp + vy * dt);
        z.push_back(z_temp + vz * dt);
        
        ax.push_back(0.0);
        ay.push_back(0.0);
        az.push_back(0.0);
        
        massa.push_back(m);
        total_mass += m;
    }
    data.close();
    
    const int N = x.size();

    cout << "N = " << N << " estrelas (ficheiro: " << input_file << ")" << endl;
    cout << "Massa total = " << total_mass << " M☉" << endl;
    
    ofstream amount(output_file);
    
    cout << "Iterações: " << iter << endl;
    cout << "dt = " << dt << " anos" << endl;
    cout << "Tempo total = " << (iter * dt / 1e6) << " milhões de anos" << endl;
    
    vector<double> x_new(N);
    vector<double> y_new(N);
    vector<double> z_new(N);
    
    for (int i = 0; i < iter; i++) {
        for (int j = 0; j < N; j++) {
            ax[j] = 0.0;
            ay[j] = 0.0;
            az[j] = 0.0;
        }
        
        for (int j = 0; j < N; j++) {
            for (int k = j+1; k < N; k++) {
                double dx = x[k] - x[j];
                double dy = y[k] - y[j];
                double dz = z[k] - z[j];
                
                double r_sq = dx*dx + dy*dy + dz*dz;
                double r = sqrt(r_sq);
                double r_inv3 = 1.0 / (r * r_sq);
                
                double fx = G * dx * r_inv3;
                double fy = G * dy * r_inv3;
                double fz = G * dz * r_inv3;
                
                ax[j] += fx * massa[k];
                ay[j] += fy * massa[k];
                az[j] += fz * massa[k];
                
                ax[k] -= fx * massa[j];
                ay[k] -= fy * massa[j];
                az[k] -= fz * massa[j];
            }
        }

        for (int j = 0; j < N; j++) {
            x_new[j] = 2.0 * x[j] - x_old[j] + ax[j] * dt_sq;
            y_new[j] = 2.0 * y[j] - y_old[j] + ay[j] * dt_sq;
            z_new[j] = 2.0 * z[j] - z_old[j] + az[j] * dt_sq;
        }
        
        double cm_x = 0.0, cm_y = 0.0, cm_z = 0.0;
        for (int j = 0; j < N; j++) {
            cm_x += x_new[j] * massa[j];
            cm_y += y_new[j] * massa[j];
            cm_z += z_new[j] * massa[j];
        }
        cm_x /= total_mass;
        cm_y /= total_mass;
        cm_z /= total_mass;
        
        int count = 0;
        for (int j = 0; j < N; j++) {
            double dx = x_new[j] - cm_x;
            double dy = y_new[j] - cm_y;
            double dz = z_new[j] - cm_z;
            
            double dist_sq = dx*dx + dy*dy + dz*dz;
            
            if (dist_sq < r_t_sq) {
                count++;
            }
        }
        
        double tempo = i * dt;
        double fracao = (double)count / N;
        amount << tempo << "\t" << fracao << "\t" << count << endl;
        
        for (int j = 0; j < N; j++) {
            x_old[j] = x[j];
            y_old[j] = y[j];
            z_old[j] = z[j];
            
            x[j] = x_new[j];
            y[j] = y_new[j];
            z[j] = z_new[j];
        }
        
        if (i % 100 == 0) {
            cout << "Iteração " << i << "/" << iter 
                 << " - Estrelas dentro de r_t: " << count << endl;
        }
    }
    
    amount.close();
    cout << "Resultados em: " << output_file << endl;
}

void ex2(int N_estrelas) {
    string filename = to_string(N_estrelas) + ".dat";
    string output_filename = "ex2_" + to_string(N_estrelas) + ".txt";
    ex2_from_file(filename, output_filename);
} 
   
void ex3(int N, string output_file){
    const double rc = 1.0;
    const double rt = 5.0;
    const double v_p = 1.0 * VELCONV;
    
    vector<double> x, y, z, vx, vy, vz, massa;
    
    for (int i = 0; i < N; i++) {
        double m;
        double choice = drand48();
        
        if (choice < 0.0008) {
            m = gerar_massa(0.01, 0.08, 0.3);
        } else if (choice < 0.005) {
            m = gerar_massa(0.08, 0.5, 1.3);
        } else {
            m = gerar_massa(0.5, 100.0, 2.3);
        }
        
        double r = gerar_raio_king(rc, rt);
        double u = drand48(), v = drand48();
        double theta = 2.0 * M_PI * u;
        double phi = acos(2.0 * v - 1.0);
        
        double x_temp = r * sin(phi) * cos(theta);
        double y_temp = r * sin(phi) * sin(theta);
        double z_temp = r * cos(phi);
        
        double vx_temp = randn() * v_p / sqrt(2.0);
        double vy_temp = randn() * v_p / sqrt(2.0);
        double vz_temp = randn() * v_p / sqrt(2.0);
        
        x.push_back(x_temp);
        y.push_back(y_temp);
        z.push_back(z_temp);
        vx.push_back(vx_temp);
        vy.push_back(vy_temp);
        vz.push_back(vz_temp);
        massa.push_back(m);
    }
    
    ofstream output(output_file);
    output << "# x(pc) y(pc) z(pc) vx(pc/ano) vy(pc/ano) vz(pc/ano) massa(M☉)" << endl;
    for (int i = 0; i < N; i++) {
        output << x[i] << " " << y[i] << " " << z[i] << " "
               << vx[i] << " " << vy[i] << " " << vz[i] << " "
               << massa[i] << endl;
    }
    
    output.close();
}

int main(){
    um_ponto_um();
    return 0;
}
