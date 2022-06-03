#include<bits/stdc++.h>
using namespace std;
using namespace std::chrono;
int main(int argc, char *argv[])
{
    int N=atoi(argv[1]);
    int T=atoi(argv[2]);
    
    /*double a[N][N];        //{{1,2,3},{3,1,4},{5,3,1}};
    double b[N][N];        //{{1,2,3},{3,1,4},{5,3,1}};
    double mult[N][N]; 
    double p[N], u[N][N], l[N][N];*/
    int *p;
    double **a,**b, **l, **u,**mult;
        a=(double **)(new double[N]);
        b=(double **)(new double[N]);
        mult=(double **)(new double[N]);
        l=(double **)(new double[N]);
        u=(double **)(new double[N]);
        p=new int[N];
        #pragma omp parallel for num_threads(T)
        for(int i=0; i<N; i++)
        {
            a[i]=new double[N];
            b[i]=new double[N];
            mult[i]=new double[N];
            l[i]=new double[N];
            u[i]=new double[N];
        }
    srand(time(0));
    #pragma omp parallel for collapse(2) num_threads(T)
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        { 
            a[i][j] =b[i][j]=int(rand()%100);
        }
    }
    #pragma omp parallel for collapse(2) num_threads(T) 
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            u[i][j] = 0.0;
            if(i==j)
                l[i][j] = 1.0;
            else
                l[i][j] = 0.0;
        }
    }
    #pragma omp parallel for num_threads(T)
    for(int i=0; i<N; i++)
    {
        p[i] = i;
    }
    auto start = high_resolution_clock::now();
    int kd;
    for(int k=0; k<N; k++)
    {
        double mx = 0.0;
        for(int i=k; i<N; i++)
        {
            if(mx < abs(a[i][k]))
            {
                mx = a[i][k];
                kd = i;
            }
        }
        if(mx == 0.0)
        {
            cout<<"Matrix is singular"<<endl;
            return 0;
        }
        swap(p[k], p[kd]);
        #pragma omp parallel for num_threads(T)
        for(int i=0;i<N;i++)
        {
            swap(a[k][i],a[kd][i]);
        }
        #pragma omp parallel for num_threads(T)
        for(int i=0; i<=k-1; i++)
        {
            swap(l[k][i], l[kd][i]);
        }

        u[k][k] = a[k][k];
        #pragma omp parallel for num_threads(T)
        for(int i=k+1; i<N; i++)
        {
            l[i][k] = a[i][k]/u[k][k];
            u[k][i] = a[k][i];
        }
        #pragma omp parallel for collapse(2) num_threads(T)
        for(int i=k+1; i<N; i++)
        {
            for(int j=k+1; j<N; j++)
            {
                a[i][j] = a[i][j] - (l[i][k]*u[k][j]);
            }
        }
    }     
    cout<<endl;
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);

    cout<<"Matrix A : "<<endl;

    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            cout<<b[i][j]<<" ";
        }
        cout<<endl;
    } 
   for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            for(int k = 0; k < N; ++k)
            {
                mult[i][j] += l[i][k] * u[k][j];
            }
    cout<<endl<<"Multiplication Matrix LU : "<<endl;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            cout<<mult[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl<<"Pivot P :"<<endl;
    for(int i=0; i<N; i++)
    {
        cout<<p[i]<<" ";
    }
    cout<<endl;
    cout<<endl<<"Matrix L : "<<endl;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            cout<<l[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl<<"Matrix U :"<<endl;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            cout<<u[i][j]<<" ";
        }
        cout<<endl;
    }
    cout <<endl<<"Time taken: " << duration.count() << " micro seconds" <<endl;
    cout<<endl;
 return 0;
}