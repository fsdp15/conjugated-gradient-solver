#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include <ctype.h>
#include <likwid.h>

#define USE_LIKWID_PERFCTR

int generateRandomDiagonal( unsigned int N, unsigned int k, unsigned int nBandas, double *diag ){
  int i;
  
  if ( !diag || N < 3 || nBandas > N/2 || k < 0 || k > nBandas )
    return (-1);

  /* garante valor dominante para diagonal principal */
  double fator = (k == 0) ? ((double)(nBandas-1)) : (0.0);

  double invRandMax = 1.0 / (double)RAND_MAX;

  for (i=0; i < N-k; ++i){
    diag[i] = fator + (double)rand() * invRandMax;;
  }

  return (0);
}

double timestamp(void){
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

void define_parametros(char *argv[], int argc, int *n, int *nBandas, int *maxIter, double *tol, FILE **arq) {
    //função para definir as variaveis e os parametros de argv
    if ((argc < 5) || (argc > 9)) {
        fprintf(stderr, "Erro: número de parâmetros incorreto.\n");
        exit(1);
    }

    int ver = 0;
    *n = atoi(argv[1]);
    *nBandas = (atoi(argv[2]) - 1) / 2;

    if (*n <= 0 || ((*nBandas == 0) && ((char)*argv[2] != '0'))) {
        fprintf(stderr, "Erro: parâmetros obrigatórios estão incorretos.\n");
        exit(1);
    }

    if(strcmp(argv[3], "-i") == 0) *maxIter = atoi(argv[4]);
    else {
        *maxIter = *n;
        ver -= 2;
    }
    if(strcmp(argv[5 + ver], "-t") == 0) sscanf(argv[6 + ver],"%lf",tol);
    else {
        *tol = (double)0;
        ver -= 2;
    }

    *arq = fopen(argv[8 + ver], "w");
    if (*arq == NULL) {
        fprintf(stderr, "Erro: não foi passado arquivo de saída.\n");
        exit(1);
    }
}

void gera_valores (double *diagonais, double *residuo, double *v, double *aux, unsigned int n, unsigned int nBandas) {
    int i;
    int percorre = n;

    if (generateRandomDiagonal(n, 0, 2*nBandas + 1, &diagonais[0]) == -1) {
        fprintf(stderr, "Erro: tamanho de banda inválido.\n");
        exit(1);
        }

    for(i=1;i<=nBandas;i++) { // Gerando matriz simétrica
        generateRandomDiagonal(n, i, 2*nBandas + 1, &diagonais[percorre]);
        percorre += n; // Percorre o vetor como se fosse uma matriz bidimensional
    }

    int aux2 = n - 4;

    for(i=0;i<aux2;i+=4) { // Residuo recebe os termos independentes, v recebe resíduo e aux = r^t (r transposto) * r
        residuo[i] = 4 * M_PI * M_PI *(sin(2 * M_PI * (i * M_PI)/n) + sin(2*M_PI*(M_PI - (i * M_PI)/n)));
        residuo[i+1] = 4 * M_PI * M_PI *(sin(2 * M_PI * ((i+1) * M_PI)/n) + sin(2*M_PI*(M_PI - ((i+1) * M_PI)/n)));
        residuo[i+2] = 4 * M_PI * M_PI *(sin(2 * M_PI * ((i+2) * M_PI)/n) + sin(2*M_PI*(M_PI - ((i+2) * M_PI)/n)));
        residuo[i+3] = 4 * M_PI * M_PI *(sin(2 * M_PI * ((i+3) * M_PI)/n) + sin(2*M_PI*(M_PI - ((i+3) * M_PI)/n)));

        v[i] = residuo[i];
        v[i+1] = residuo[i+1];
        v[i+2] = residuo[i+2];
        v[i+3] = residuo[i+3];

        *aux += residuo[i] * residuo[i];
        *aux += residuo[i+1] * residuo[i+1];
        *aux += residuo[i+2] * residuo[i+2];
        *aux += residuo[i+3] * residuo[i+3];
    }

    for (; i < n; i++) {
        residuo[i] = 4 * M_PI * M_PI *(sin(2 * M_PI * (i * M_PI)/n) + sin(2*M_PI*(M_PI - (i * M_PI)/n)));
        v[i] = residuo[i];
        *aux += residuo[i] * residuo[i];    
    }
}

void multiplica_matriz_vetor(double *destino, double *A, double *vetor, unsigned int n, unsigned int nBandas) {
  int i, j, aux, aux2, aux3;
  aux3 = n - 4;

    likwid_markerStartRegion("mult_matriz");
    for (i = 0; i < aux3; i += 4) {
        destino[i] = A[i] * vetor[i];
        destino[i+1] = A[i+1] * vetor[i+1];
        destino[i+2] = A[i+2] * vetor[i+2];
        destino[i+3] = A[i+3] * vetor[i+3];
    }

    for(; i < n; i++) {
        destino[i] = A[i] * vetor[i];
    }

    for (i = 1; i <= nBandas; i++) {
        aux = i;
        aux2 = n-i-1;
        aux3 = i*n;
        for (j = 0; j < n-i-4; j += 4) {
            destino[aux] += A[aux3 + j] * vetor[j];
            destino[aux+1] += A[aux3 + j + 1] * vetor[j + 1];
            destino[aux+2] += A[aux3 + j + 2] * vetor[j + 2];
            destino[aux+3] += A[aux3 + j + 3] * vetor[j + 3];
            aux += 4;
            destino[aux2] += A[aux3 + aux2] * vetor[aux2+i];
            destino[aux2-1] += A[aux3 + aux2 - 1] * vetor[aux2+i - 1];
            destino[aux2-2] += A[aux3 + aux2 - 2] * vetor[aux2+i - 2];
            destino[aux2-3] += A[aux3 + aux2 - 3] * vetor[aux2+i - 3];
            aux2 -= 4;
        }

        for (; j < n-i; j++) {
            destino[aux] += A[aux3 + j] * vetor[j];
            aux++;
            destino[aux2] += A[aux3 + aux2] * vetor[aux2+i];
            aux2--;
        }
    } 
    likwid_markerStopRegion("mult_matriz");
}

void calculaS(double aux, double *v, double *z, double *minimiza, unsigned int n) { // S = minimiza
    double mult; // Variável para armazenar o resultado de v^t * z
    int i, aux2;
    aux2 = n - 4;

    mult = (double)0;

    likwid_markerStartRegion("mult_vetor");  

    for (i = 0; i < aux2; i += 4) {
        mult += v[i] * z[i];
        mult += v[i+1] * z[i+1];
        mult += v[i+2] * z[i+2];
        mult += v[i+3] * z[i+3];
    }

    for (; i < n; i++) 
        mult += v[i] * z[i];

    likwid_markerStopRegion("mult_vetor");

    if (mult != (double)0)
        *minimiza = aux / mult;
    else {
        fprintf(stderr, "Erro: divisão por zero. Diminua o número de iterações ou aumente a tolerância.\n");
        exit(1);
    }
}

void calcula_residuo(double *residuo, double *z, double *parcial, double *v, double minimiza, unsigned int n) {
    int i;

    for(i=0;i<n;i++) {
        residuo[i] -= minimiza * z[i];
        parcial[i] += minimiza * v[i];
        if (parcial[i] != parcial[i]) {
            fprintf(stderr, "Erro: overflow (NaN). Diminua o número de iterações ou aumente a tolerância.\n");
            exit(1);
      }
    }
}

void define_outros(double *aux, double *v, double *residuo, unsigned int n) {
    int i, aux3;
    double m;
    double aux1 = (double)0;
    aux3 = n-4;

    for(i=0;i<aux3;i+=4) {
        aux1 += residuo[i] * residuo[i];  
        aux1 += residuo[i+1] * residuo[i+1];  
        aux1 += residuo[i+2] * residuo[i+2];  
        aux1 += residuo[i+3] * residuo[i+3];  
    }

    for (; i < n; i++) 
        aux1 += residuo[i] * residuo[i];  

    if (*aux != (double)0)
        m = aux1 / *aux;
    else {
        fprintf(stderr, "Erro: divisão por zero. Diminua o número de iterações ou aumente a tolerância.\n");
        exit(1);
    }

    *aux = aux1;

    for (i=0;i<aux3;i+=4) {
        v[i] = residuo[i] + (m * v[i]);
        v[i+1] = residuo[i+1] + (m * v[i+1]);
        v[i+2] = residuo[i+2] + (m * v[i+2]);
        v[i+3] = residuo[i+3] + (m * v[i+3]);
    }

    for (; i < n; i++) 
        v[i] = residuo[i] + (m * v[i]);

}


void imprime_resultados(double tempoMetodoMin, double tempoMetodoMedia, double tempoMetodoMax, double tempoResiduoMin, double tempoResiduoMedia,
                        double tempoResiduoMax, double *norma, double *parcial, FILE* arq, unsigned int n, int p) {
    //funcao que imprime os resultados do programa em um arquivo de saida
    int i;

    fprintf(arq,"###########\n");
    fprintf(arq,"# Tempo Método CG: %g %g %g\n", tempoMetodoMin, tempoMetodoMedia, tempoMetodoMax);
    fprintf(arq,"# Tempo Resíduo: %g %g %g\n#\n", tempoResiduoMin, tempoResiduoMedia, tempoResiduoMax);
    fprintf(arq,"Norma Euclidiana do Residuo e Erro aproximado\n");
    fprintf(arq,"# i=%d: %g %g\n", 1, norma[0], (double)0);
    for(i=1;i<p;i++)
        fprintf(arq,"# i=%d: %g %g\n", i+1, norma[i], fabs(norma[i] - norma[i-1]));
    fprintf(arq,"###########\n%d\n", n);
    for (i = 0; i < n; i++)
        fprintf(arq, "%.14g ", parcial[i]);
    free(norma);
    free(parcial);
    free(arq);
}


int main(int argc, char *argv[]){
    likwid_markerInit();
    int i, p, maxIter;
    int n, nBandas; // Variáveis para armanezar os parâmetros digitados no terminal
    double tol;
    // variaveis para cálculos de tempo
    double tempoMetodoMedia, tempoResiduoMedia = (double)0,
           tempoM, tempoR;
    double tempoResiduoMax = (double)0, tempoResiduoMin = (double)10,
           tempoMetodoMax = (double)0, tempoMetodoMin = (double)10;

    FILE *arq = NULL;
    define_parametros(argv, argc, &n, &nBandas, &maxIter, &tol, &arq);

    double *diagonais; //Armaneza as diagonais de forma linear em um vetor unidimensional
    double *parcial = calloc(n, sizeof(double)); // Vetor Xk para armanezar o resultado parcial xnBandas de cada iteração
    double *z, *v; // Vetores auxiliares
    double *residuo, *norma, erro; // Vetores para armazenar o resíduo, norma e erro
    // Não é utilizado um vetor para armazenar os termos independentes pois no inicio da execução o vetor residuo ja recebe os termos independentes
    double minimiza = (double)0; // Valor s que minimiza a função
    double aux = (double)0; 

    srand(20162);
    diagonais = malloc(sizeof(double) * (nBandas + 1) * n); // Alocando vetores
    z = malloc(sizeof(double) * n);
    v = malloc(sizeof(double) * n);
    parcial = malloc(sizeof(double) * n);
    residuo = malloc(sizeof(double) * n);
    norma = malloc(sizeof(double) * n);

    gera_valores(diagonais, residuo, v, &aux, n, nBandas); // Gerando a matriz de banda
    erro = tol + 1.0;
    tempoMetodoMedia = timestamp();

    for (p = 0; ((p < maxIter) && (erro > tol)); p++) { //inicio da execução do metodo
        tempoM = timestamp();
        tempoR = timestamp();
        multiplica_matriz_vetor(z, diagonais, v, n, nBandas); // z = Av
        calculaS (aux, v, z, &minimiza, n); // Calcula a variável s (minimiza) que é usada para minimizar a solução do sistema
        calcula_residuo(residuo, z, parcial, v, minimiza, n); // residuo r = minimiza * z

        tempoR = timestamp() - tempoR;

        if((tempoR) < tempoResiduoMin)
            tempoResiduoMin = tempoR;
        if((tempoR) > tempoResiduoMax)
            tempoResiduoMax = tempoR;
        tempoResiduoMedia += tempoR;
        norma[p] = (double)0;
        for(i=0;i<n;i++)
            norma[p] += (residuo[i]*residuo[i]);
        norma[p] = sqrt(norma[p]); // Calculando norma euclidiana do resíduo
        if(p!=0)
            erro = fabs(norma[p] - norma[p-1]); // Calculando erro aproximado absoluto

        define_outros(&aux, v, residuo, n); // Atualiza vetor auxiliar v e variável aux para a próxima iteração

        tempoM = timestamp() - tempoM;
        if((tempoM) < tempoMetodoMin)
            tempoMetodoMin = tempoM;
        if((tempoM) > tempoMetodoMax)
            tempoMetodoMax = tempoM;

    } //fim da execução do método
    likwid_markerClose();
    free(diagonais);
    free(z);
    free(v);
    free(residuo);

    tempoMetodoMedia = timestamp() - tempoMetodoMedia ; //calcula tempo que levou para executar o metodo

    tempoResiduoMedia /= p; //calcula tempo medio de execução do residuo
    tempoMetodoMedia /= p;  //calculo do tempo medio de execução do método
    imprime_resultados(tempoMetodoMin, tempoMetodoMedia, tempoMetodoMax, tempoResiduoMin,
                       tempoResiduoMedia, tempoResiduoMax, norma, parcial, arq, n, p);
    return 0;

}
