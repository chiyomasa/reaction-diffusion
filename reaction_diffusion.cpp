#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define X 800
#define Y 800
#define h 0.1
#define roop 10000
#define dt 0.001
#define DU 0.035
#define DV 1.0
#define U0 1.0
#define V0 14.5
#define W0 0.0

typedef int unit;
typedef double UNIT;

typedef struct
{
	UNIT NUM;
	UNIT NUMR;
	UNIT NUMD;
	UNIT SUM;
}Calculation;

UNIT REACT_U(UNIT u, UNIT v);
UNIT REACT_V(UNIT u, UNIT v);
UNIT INACTIVE(UNIT u, UNIT v);
UNIT Random();
unit P_N();

int main()
{
	unit x, y, t;
	unit count, output_uw, output;
	unit pn;
	UNIT CV;
	Calculation u[X + 1][Y + 1], v[X + 1][Y + 1], w[X + 1][Y + 1];
	UNIT CU[X + 1][Y + 1];


	FILE *fp,*gp;
	char FILEname[50];

	count = 0;
	output_uw = roop / 200;
	output = roop/10;

	CV = DV*dt / (h*h);

	//初期値の設定

	for (y = 0; y <= Y; y++)
	{
		for (x = 0; x <= X; x++)
		{
			if ((pow(x - 400, 2) + pow(y - 400, 2)) <= pow(15, 2))
			{
				u[x][y].NUM = U0;
				v[x][y].NUM = V0;
				w[x][y].NUM = W0;
			}
			else
			{
				u[x][y].NUM = 0.0;
				v[x][y].NUM = V0;
				w[x][y].NUM = W0;
			}
		}
	}

	//初期値のファイル出力

	sprintf_s(FILEname, "initial.dat");
	fopen_s(&fp, FILEname, "w");
	for (y = 0; y <= Y; y++)
	{
		for (x = 0; x <= X; x++)
		{
			fprintf_s(fp, "%d %d %lf\n", x, y, u[x][y].NUM + w[x][y].NUM);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	//反応拡散方程式

	//拡散係数（バクテリア）の設定

	srand((unsigned int)time(NULL));
	for (y = 0; y <= Y; y++)
	{
		for (x = 0; x <= X; x++)
		{
			pn = P_N();
			if (pn >= 50)
			{
				CU[x][y] = DU + 0.01*Random();
			}
			else
			{
				CU[x][y] = DU - 0.01*Random();
			}

		}
	}

	//拡散係数の出力

	sprintf_s(FILEname, "coefficientU.dat");
	fopen_s(&fp, FILEname, "w");

	for (y = 0; y <= Y; y++)
	{
		for (x = 0; x <= X; x++)
		{
			fprintf_s(fp, "%d %d %lf\n", x, y, CU[x][y]);
		}
		fprintf_s(fp,"\n");
	}

	fclose(fp);

	for (t = 0; t <= roop; t++)
	{

		//境界条件（ノイマン）

		for (x = 0; x <= X; x++)
		{
			u[x][0].NUM = u[x][1].NUM;
			v[x][0].NUM = v[x][1].NUM;
			w[x][0].NUM = w[x][1].NUM;

			u[x][Y].NUM = u[x][Y - 1].NUM;
			v[x][Y].NUM = v[x][Y - 1].NUM;
			w[x][Y].NUM = w[x][Y - 1].NUM;
		}

		for (y = 0; y <= Y; y++)
		{
			u[0][y].NUM = u[1][y].NUM;
			v[0][y].NUM = v[1][y].NUM;
			w[0][y].NUM = w[1][y].NUM;

			u[X][y].NUM = u[X - 1][y].NUM;
			v[X][y].NUM = v[X - 1][y].NUM;
			w[X][y].NUM = w[X - 1][y].NUM;
		}

		//反応項の計算

		for (y = 1; y <= Y - 1; y++)
		{
			for (x = 1; x <= X - 1; x++)
			{
				u[x][y].NUMR = REACT_U(u[x][y].NUM, v[x][y].NUM);
				v[x][y].NUMR = REACT_V(u[x][y].NUM, v[x][y].NUM);
			}
		}

		//拡散項の計算
		//x=1~X-1,y=1~Y-1
		for (y = 1; y <= Y - 1; y++)
		{
			for (x = 1; x <= X - 1; x++)
			{
				u[x][y].NUMD = u[x][y].NUM + CU[x][y] * dt / (h*h)*(u[x + 1][y].NUM + u[x - 1][y].NUM + u[x][y + 1].NUM + u[x][y - 1].NUM - 4.0*u[x][y].NUM);
				v[x][y].NUMD = v[x][y].NUM + CV*(v[x + 1][y].NUM + v[x - 1][y].NUM + v[x][y + 1].NUM + v[x][y - 1].NUM - 4.0*v[x][y].NUM);
			}
		}

		//反応拡散方程式

		for (y = 1; y <= Y - 1; y++)
		{
			for (x = 1; x <= X - 1; x++)
			{
				u[x][y].SUM = u[x][y].NUMR + u[x][y].NUMD - INACTIVE(u[x][y].NUM, v[x][y].NUM);
				v[x][y].SUM = v[x][y].NUMR + v[x][y].NUMD;
				w[x][y].SUM = w[x][y].NUM + INACTIVE(u[x][y].NUM, v[x][y].NUM);
			}
		}

		//値の更新

		for (y = 1; y <= Y - 1; y++)
		{
			for (x = 1; x <= X - 1; x++)
			{
				u[x][y].NUM = u[x][y].SUM;
				v[x][y].NUM = v[x][y].SUM;
				w[x][y].NUM = w[x][y].SUM;
			}
		}

		//ファイルへの出力

		if (t%output == 0)
		{

			//u
			sprintf_s(FILEname, "colony[%d]_u.dat", t);
			fopen_s(&fp, FILEname, "w");
			for (y = 0; y <= Y; y++)
			{
				for (x = 0; x <= X; x++)
				{
					fprintf_s(fp, "%d %d %lf\n", x, y, u[x][y].NUM);
				}
				fprintf_s(fp, "\n");
			}
			fclose(fp);

			gp = _popen("pgnuplot.exe", "w");
			fprintf_s(gp, "set terminal png\n");
			fprintf_s(gp, "set size square\n");
			fprintf_s(gp, "set view map\n");
			fprintf_s(gp, "set palette cubehelix start 1 cycles 0 saturation 2 negative\n");
			fprintf_s(gp, "set output \"colony[%d]_u.png\"\n", t);
			fprintf_s(gp, "splot \"colony[%d]_u.dat\" with pm3d\n", t);
			fflush(gp);
			_pclose(gp);


			//v
			sprintf_s(FILEname, "colony[%d]_v.dat", t);
			fopen_s(&fp, FILEname, "w");
			for (y = 0; y <= Y; y++)
			{
				for (x = 0; x <= X; x++)
				{
					fprintf_s(fp, "%d %d %lf\n", x, y, v[x][y].NUM);
				}
				fprintf_s(fp, "\n");
			}
			fclose(fp);

			gp = _popen("pgnuplot.exe", "w");
			fprintf_s(gp, "set terminal png\n");
			fprintf_s(gp, "set size square\n");
			fprintf_s(gp, "set view map\n");
			fprintf_s(gp, "set palette cubehelix start 1 cycles 0 saturation 2 negative\n");
			fprintf_s(gp, "set output \"colony[%d]_v.png\"\n", t);
			fprintf_s(gp, "splot \"colony[%d]_v.dat\" with pm3d\n", t);
			fflush(gp);
			_pclose(gp);

			//w
			sprintf_s(FILEname, "colony[%d]_w.dat", t);
			fopen_s(&fp, FILEname, "w");
			for (y = 0; y <= Y; y++)
			{
				for (x = 0; x <= X; x++)
				{
					fprintf_s(fp, "%d %d %lf\n", x, y, w[x][y].NUM);
				}
				fprintf_s(fp, "\n");
			}
			fclose(fp);

			gp = _popen("pgnuplot.exe", "w");
			fprintf_s(gp, "set terminal png\n");
			fprintf_s(gp, "set size square\n");
			fprintf_s(gp, "set view map\n");
			fprintf_s(gp, "set palette cubehelix start 1 cycles 0 saturation 2 negative\n");
			fprintf_s(gp, "set output \"colony[%d]_w.png\"\n", t);
			fprintf_s(gp, "splot \"colony[%d]_w.dat\" with pm3d\n", t);
			fflush(gp);
			_pclose(gp);
		}
		if ((t%output_uw) == 0)
		{
			sprintf_s(FILEname, "colony[%d]_uw.dat", t);
			fopen_s(&fp, FILEname, "w");
			for (y = 0; y <= Y; y++)
			{
				for (x = 0; x <= X; x++)
				{
					fprintf_s(fp, "%d %d %lf\n", x, y, u[x][y].NUM + w[x][y].NUM);
				}
				fprintf_s(fp, "\n");
			}
			fclose(fp);

			gp = _popen("pgnuplot.exe", "w");
			fprintf_s(gp, "set terminal png\n");
			fprintf_s(gp, "set size square\n");
			fprintf_s(gp, "set view map\n");
			fprintf_s(gp, "set palette cubehelix start 1 cycles 0 saturation 2 negative\n");
			fprintf_s(gp, "set output \"colony[%d]_uw.png\"\n", t);
			fprintf_s(gp, "splot \"colony[%d]_uw.dat\" with pm3d\n", t);
			fflush(gp);

			_pclose(gp);
		}

		count++;

		if (count % 1000 == 0)
		{
			printf("%d回目\n", count);
		}

	}

	return 0;

}
UNIT REACT_U(UNIT u, UNIT v)
{
	return 20.0*dt*u*v;
}

UNIT REACT_V(UNIT u, UNIT v)
{
	return -dt*u*v;
}

UNIT INACTIVE(UNIT u, UNIT v)
{
	return 2400 * dt*u / ((1.0 + u)*(1.0 + v));
}

UNIT Random()
{
	return (double)rand() / RAND_MAX;
}

unit P_N()
{
	return rand() % 100 + 1;
}
