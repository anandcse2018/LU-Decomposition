#include<bits/stdc++.h>
using namespace std;
using namespace std::chrono;

pthread_mutex_t mutex1=PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t barrier;

int counter=0;
int T=4;
int N=4;
struct matrix
{
	int *p;
	double **a, **l, **u, **b;
	matrix()
	{	
		a=(double **)(new double[N]);
		b=(double **)(new double[N]);
		l=(double **)(new double[N]);
		u=(double **)(new double[N]);
		p=new int[N];

		for(int i=0; i<N; i++)
		{
			a[i]=new double[N];
			b[i]=new double[N];
			l[i]=new double[N];
			u[i]=new double[N];
			l[i][i] = 1;
		}
	}
};
struct algo
{
	matrix *datap;
	int k;
	int kd;
	algo(matrix *data, int z, int y)
	{
		datap = data;
		k = z;
		kd = y;
	}
};
void *in(void *arg)
{
	matrix *datap = (matrix *)arg;
	double **a = datap->a;
	double **b = datap->b;
	int* p = datap->p;
	srand(time(NULL));
	for(int i=0; i<N; i++)
    {
    	p[i] = i;
        for(int j=0; j<N; j++)
        { 
            a[i][j]=b[i][j]=int(rand()%100);
        }
    }
}
void *process(void *arg)
{
	algo *datap_p = (algo *)arg;
	double **a = datap_p->datap->a;
	double **l = datap_p->datap->l;
	double **u = datap_p->datap->u;
	int k = datap_p->k;
	int kd = datap_p->kd;

	pthread_mutex_lock( &mutex1 );
	int tid = counter++;
	pthread_mutex_unlock( &mutex1 );
	
	int c1=ceil(((k*(1.0))/T));
	int c2=ceil((((N-k-1)*(1.0))/T));
	int start1=tid*c1;
	int last1=start1+c1;
	last1=last1>k?k:last1;
	int start2=(k+1)+tid*c2;
	int last2=start2 + c2;
	last2=last2>N?N:last2;
	
	for(int i=start1; i<last1; i++)
	{
		double temp_l = l[k][i];
		l[k][i] = l[kd][i];
		l[kd][i] = temp_l;
	}
	double u_k_k = u[k][k];
	for(int i=start2; i<last2; i++)
	{
		l[i][k] = a[i][k]/u_k_k;
		u[k][i] = a[k][i];
	}	
	pthread_barrier_wait(&barrier);
	for(int i=start2; i<last2; i++)
	{
		double operand1 = l[i][k];
		for(int j=k+1; j<N; j++)
			a[i][j] -= operand1*u[k][j];
	}
}
int get_max(double** a, int c, int r1, int r2)
{
	double max=0;
	int index=-1;
	for(int i= r1; i<= r2; i++)
	{
		double temp = abs(a[i][c]);
		if(temp > max)
		{
			max = temp;
			index = i;
		}
	}
	if(max == 0)
	{
		cout<<"singular Matrix"<<endl;
		return 0;
	}
	return index;
}
void print_matrix(double **a, int size)
{
	for(int i=0; i<size; i++)
	{
		for(int j=0; j<size; j++)	cout << a[i][j] << " ";
		cout << endl;
	}
}
void print_vector(int *vec, int size)
{
	for(int i=0; i<size; i++)	cout << vec[i] << " ";
	cout << endl;
}
void multiply(double** u, double** l, int size)
{
	for(int i=0; i<size; i++)
	{
		for(int j=0; j<size; j++)
		{
			double temp = 0;
			for(int k=0; k<size; k++)
			{  
				temp += l[i][k]*u[k][j];
			}
			cout << temp << " ";
		}
		cout << endl;
	}
}
int main(int argc, char *argv[])
{
	pthread_barrier_init(&barrier, NULL, T);
	sscanf(argv[1], "%d", &N);
	sscanf(argv[2], "%d", &T);	

	matrix data = matrix();
	matrix *datap = &data;	
	pthread_barrier_init(&barrier, NULL, T);
	pthread_t threads[T];

	for(int i=0; i<T; i++)
		pthread_create(&threads[i], NULL, in, (void *)datap);

	for(int i=0; i<T; i++)
		pthread_join(threads[i], NULL);
    auto start1 = high_resolution_clock::now();
	for(int k=0; k<N; k++)
	{
		int kd = get_max(data.a, k, k, N-1);
		if(kd == -1)
		{
			cout << "Singular Matrix\n";
			return 0;
		}
		int temp_p = data.p[k];
		data.p[k] = data.p[kd];
		data.p[kd] = temp_p;
		
		double* temp_row_pointer = data.a[k];
		data.a[k] = data.a[kd];
		data.a[kd] = temp_row_pointer;
		
		data.u[k][k] = data.a[k][k];
		
		algo data_p = algo(&data, k, kd);
		algo *datap_p = &data_p;
		for(int i=0; i<T; i++)
			pthread_create(&threads[i], NULL, process, (void *)datap_p);

		for(int i=0; i<T; i++)
			pthread_join(threads[i], NULL);
		counter = 0;
	}	
	auto end = high_resolution_clock::now();	
	auto duration_time = duration_cast<microseconds>(end - start1);

	cout<<endl<<"Original Matrix :"<<endl;
    print_matrix(data.b,N);
    cout<<endl<<"Multiplication Matrix LU :"<<endl;
    multiply(data.u,data.l,N);
    cout<<endl<<"Permutation Matrix :"<<endl;
    print_vector(data.p, N);
    cout<<endl<<"Lower Triangular Matrix :"<<endl;
	print_matrix(data.l,N);
	cout<<endl<<"Upper Triangular Matrix :"<<endl;
    print_matrix(data.u,N);
    cout<<endl<<"Time taken: "<<duration_time.count()<<" microseconds" <<endl;
    cout<<endl;
	return 0;
}