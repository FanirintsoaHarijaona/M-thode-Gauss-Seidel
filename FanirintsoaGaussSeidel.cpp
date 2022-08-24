#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;

void initialiserMatrice(vector<vector<float>> &matrice,vector<float> &second,int &dim);
vector<float> resolutionGaussSeidel(vector<vector<float>>& A,vector<float>&B,int &dim);
float sommeElementsMatrice(vector<float>matrice);
void afficherMatrice(vector<vector<float>>& matrice,int& dim);
void afficherVecteur(vector<float>&matrice,int  &dim);

int main(){
    vector<vector<float>> A;
    vector<float>B;
    vector<float> solution;
    int dim = 0;
    cout<<"Résolvons le système d'équations ci-dessous à l'aide de la méthode de Gauss-Seidel:"<<endl;
    initialiserMatrice(A,B,dim);
    solution = resolutionGaussSeidel(A,B,dim);
    cout<<"Le système peut se traduire sous forme matricielle en une équation A.X = B\nLa matrice au premier membre de l'équation est A:"<<endl;
    afficherMatrice(A,dim);
    cout<<"\nLe vecteur second membre de l'équation est B:"<<endl;
    afficherVecteur(B,dim);
    cout<<"L'équation A.X = B avec pour solution X"<<endl;
    afficherVecteur(solution,dim);
    return 0;
}

void afficherVecteur(vector<float>&matrice,int  &dim){
    for(int i =0 ; i<dim ; i++){
        cout<<'\t'<<matrice[i]<<endl;
    }
}

float sommeElementsMatrice(vector<float>matrice){
    float result(0);
    int dim = matrice.size();
    for(int i =0 ; i<dim ; i++){
        result += fabs(matrice[i]);
    }
    return result;
}

vector<float> resolutionGaussSeidel(vector<vector<float>>& A,vector<float>&B, int &dim){
    vector <float> x, residus;
//variable pour compter l'itération
    int iteration(1);
    float epsilone = pow(10,-6),sommeRes(1);
//Choix de l'approximation initiale x0 
    for(int i = 0 ;i<dim;i++){
        if(i==0)
            x.push_back(1);
        x.push_back(0);
        residus.push_back(0);
    }
//calcul de la solution avec la méthode de Gauss-Seidel
    while(sommeRes>=epsilone){
    for (int i = 0; i<dim ; i++){
        float somme(0);
        for(int j =0; j<dim;j++){
            somme += A[i][j]*x[j];
        }
        residus[i] = B[i] - somme;
        x[i] = x[i] + 1/A[i][i]*residus[i];
    }
    sommeRes = sommeElementsMatrice(residus);
    iteration++;
    }
    return x;
    
}

void afficherMatrice(vector<vector<float>>& matrice,int& dim){
    for (int i=0;i<dim;i++){
       for (int j=0; j<dim;j++){
            cout<<matrice[i][j]<<"\t";
        }
       cout<<endl;
    }
}

void initialiserMatrice(vector<vector<float>> &matrice,vector<float> &second,int &dim){
//initialisation d'une variable afin d'ouvrir le fichier
    ifstream fichier("data.txt");
    if(fichier){
//la première ligne du fichier est déstinéé à être la dimension de la matrice
        fichier>>dim;
//Tant que la variable i sera inférieure à la dimension insérée ci-dessus
//nous considérerons les lignes en dessous de la première comme étant les
//coéfficients de la matrice
        for(int i=0;i<dim;i++){
            vector<float> ligne;
            for(int j=0;j<dim;j++){
                float temp(0);
                fichier>>temp;
                ligne.push_back(temp);
            }
            matrice.push_back(ligne);
        }
//le reste du fichier sera la matrice de B du système d'équation AX = B où A est 
//la matrice à triangulariser
        for(int i=0;i<dim;i++){
            float temp = 0;
            fichier>>temp;
            second.push_back(temp);
        }
        for (int i= 0;i<dim ; i++){
            cout<<"|\t";
            for(int j =0;j<dim; j++){
                cout<<matrice[i][j]<<"x["<<j<<"] +";
                if(j == dim -1)
                    cout<<matrice[i][j]<<"x["<<j<<"] =";
            }
            cout<< second[i]<<endl;
        }
    }
//En cas d'erreur d'ouverture du fichier, on affiche un message d'erreur
else{
        cout<<"Erreur lors de chargement du fichier"<<endl;
    }
}