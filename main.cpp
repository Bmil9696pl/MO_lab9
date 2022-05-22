#include <iostream>
#include "cmath"
#include "fstream"
using namespace std;

// wspolczynniki zadanego rownania rozniczkowego

const double p = 1.0, q = 0.0, r = -4.0;

// z warunow brzgowych
const double alfa = 0.0, beta = 1.0, gamma = -1.0;
const double fi = 0.0, psi = 1.0, teta = 0.0;

// przedzial
const double poczatekPrzedzialu = 0.0, koniecPrzedzialu = 1.0;

double rownanieAnalityczne(double x){
    return ( exp( 2.0-2.0*x ) - 4*exp( 4.0 - x * 2.0  )+ 4.0*exp( x * 2.0 ) - exp( 2.0 + 2.0*x ) - x + x*exp(4.0)) / ( 4.0 - 4*exp(4.0) );
}

void thomas(double *l, double *d, double *u, double *b, double *x, int N){
    double *bKopia = new double[N];
    double *dKopia = new double [N];

    dKopia[0] = d[0];
    bKopia[0] = b[0];

    //wyliczenie ni
    for (int i = 1; i < N; ++i) {
        dKopia[i] = d[i] - l[i - 1] * (u[i - 1] / dKopia[i - 1]);
    }

    //wyliczenie r
    for (int i = 1; i < N; ++i) {
        bKopia[i] = b[i] - l[i - 1] * bKopia[i - 1] / dKopia[i- 1];
    }

    x[N - 1] = bKopia[N - 1] / dKopia[N - 1];

    for (int i = N - 2; i >= 0 ; --i) {
        x[i] = (bKopia[i] - u[i] * x[i + 1]) / dKopia[i];
    }
    delete[] bKopia;
    delete[] dKopia;
}

int max(double* blad, int N){
    double max = fabs(blad[0]);
    int naj = 0;
    for (int i = 0; i < N; ++i) {
        if(fabs(blad[i]) > max){
            max = fabs(blad[i]);
            naj = i;
        }
    }
    return naj;
}

double Numerow(double h, int N) {
    double *l, *d, *u, *b, *x, *blad, xp1 = poczatekPrzedzialu, xp2 = poczatekPrzedzialu; //deklarujemy zmienne
    fstream numerow, analitycznie;
    numerow.open( "numerow.txt", ios_base::app ); //otwieramy pliki do zapisu wyników wyliczeń analitycznych oraz z Numerowa
    analitycznie.open( "analitycznie.txt", ios_base::app);
    analitycznie << scientific;
    numerow << scientific;
    cout.precision(10);
    l = new double[N]; //tworzymy wektory
    d = new double[N];
    u = new double[N];
    b = new double[N];
    x = new double[N];
    blad = new double[N];

    u[0] = alfa / h; //zgodnie z wzorami uzupełniamy wartości pierwszych wyrazów, które są inne niż pozostałe
    d[0] = beta - alfa / h;
    b[0] = -gamma;

    for(int i = 1; i < N - 1; i++) //w pętli uzupełniamy wyrazy środkowe
    {
        l[i - 1] = p / (h * h) + r / 12.0;
        d[i] = (-2.0 * p) / (h * h) + r * (10.0 / 12.0);
        u[i] = p / (h * h) + r / 12.0;
        b[i] = (xp1 + i * h - h) / 12.0 + (10.0 / 12.0) *  (xp1 + i * h) + (xp1 + i * h + h) / 12.0;
    }

    l[N - 2] = -fi / h; //końcowe wyrazy mają także specjalne dane
    d[N - 1] = -fi / h + psi;
    b[N - 1] = -teta;

    thomas( l, d, u, b, x, N ); //wykonujemy algorytm Thomasa na wyliczonych danych

    for ( int i = 0; i < N; i++ ) //obliczamy błąd bezwzględny między wyliczeniami Numerowa, a rozwiązaniem analitycznym
    {
           blad[i] = fabs( x[i] - rownanieAnalityczne( xp2 ) );
        xp2 += h;
    }
    int naj = max(blad, N); //znajdujemy największy błąd

    if(N==402) //losowa liczba
    {
        for ( int i = 0; i < N; i++ )
        {
            numerow<<xp1<<"\t"<<x[i]<<"\t"<<endl;
            analitycznie<<xp1<<"\t"<<rownanieAnalityczne( xp1)<<endl;
            xp1 += h;
        }
    }
//	cout <<h<<"\t"<<rownanieAnalityczne( xp1+naj*h )<<"\t"; //wyświetlamy wyniki
//	cout<< x[naj]<<"\t"<<blad[naj]<<endl;

    delete[] l; //usuwamy niepotrzebne wektory
    delete[] d;
    delete[] u;
    delete[] x;
    delete[] b;
    analitycznie.close(); //zamykamy pliki
    numerow.close();
    return blad[naj]; //zwracamy największy błąd
}

double konwencjonalnaTrzypunktowa(double h, int N){
    double *l, *d, *u, *b, *x, *blad, xp1 = poczatekPrzedzialu, xp2 = poczatekPrzedzialu;
    fstream konwencjonalnie;
    konwencjonalnie.open("konwencjonalnie.txt", ios_base::app);
    konwencjonalnie << scientific;
    cout.precision(10);

    l = new double[N];
    d = new double[N];
    u = new double[N];
    b = new double[N];
    x = new double[N];
    blad = new double[N];

    u[0] = alfa / h;
    d[0] = beta - alfa / h;
    b[0] = -gamma;

    for (int i = 1; i < N - 1; ++i) {
        l[i - 1] = p / (h * h) - q / (2.0 * h);
        d[i] = (-2.0 * p) / (h * h) + r;
        u[i] = p / (h * h) + q / (2.0 * h);
        b[i] = (xp1 + i * h);
    }

    l[N - 2] = -fi / h;
    d[N - 1] = -fi / h + psi;
    b[N - 1] = -teta;

    thomas(l,d,u,b,x,N);

    for(int i = 0; i < N; i++){
        blad[i] = fabs(x[i] - rownanieAnalityczne(xp2));
        xp2 += h;
    }
    int naj = max(blad, N);
    if(N == 402){
        for(int i = 0; i < N; i++){
            konwencjonalnie<<xp1<<"\t"<<x[i]<<endl;
            xp1 += h;
        }
    }
    delete[] l;
    delete[] d;
    delete[] u;
    delete[] x;
    delete[] b;
    konwencjonalnie.close();
    return blad[naj];
}

int main() {
    double h, bladNum, bladKonw;
    int N;
    fstream bledy, rzedy;
    bledy.open("bledy.txt", fstream::out);
    bledy << scientific;
    cout.precision(10);

    for(N = 2; N < 30000; N += 100){
        h = (koniecPrzedzialu - poczatekPrzedzialu) / (N-1);
        bladKonw = konwencjonalnaTrzypunktowa(h, N);
        bladNum = Numerow(h, N);
        bledy<<log10(h)<<"\t"<<log10(bladKonw)<<"\t"<< log10(bladNum)<<endl;
    }
    bledy.close();
    return 0;
}
