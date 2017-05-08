double max(double a, double b)						//����һ��max
{
	if (a > b) return a;
	return b;
}
double min(double a, double b)
{
	if (a > b) return  b;
	return a;
}
const double Pi = acos(-1.0);						//Pi�Ǧ�
const double eps = 1e-8;							//�ӽ�0��һ���������Χ
int cmp(double x)									//���㼸�������������
{
	if (fabs(x) < eps) return 0;
	if (x>0) return 1;
	return -1;
}
inline double sqr(double x) { return x*x; }			//ƽ��
typedef struct  Point			//������
{
	double x, y;
	Point(){}
	Point(double a, double b) :x(a), y(b){}
	void input(){ scanf("%lf%lf", &x, &y); }	//ûʲô�õĶ���
	/*������/��ļӼ������������˳���������ˣ�����/��ıȽϽ���������*/
	friend Point operator + (const Point &a, const Point &b){ return Point(a.x + b.x, a.y + b.y); }
	friend Point operator - (const Point &a, const Point &b){ return Point(a.x - b.x, a.y - b.y); }
	friend bool operator == (const Point &a, const Point &b){ return (!cmp(a.x - b.x)) && (!cmp(a.y - b.y)); }
	friend Point operator * (const Point &a, const double &b){ return Point(a.x * b, a.y * b); }
	friend Point operator * (const double &a, const Point &b){ return Point(a * b.x, a * b.y); }
	friend double operator * (const Point &a, const Point &b){ return (a.x * b.x + a.y * b.y); }
	friend Point operator / (const Point &a, const double &b){ return Point(a.x / b, a.y / b); }
	double norm(){ return sqrt(sqr(x) + sqr(y)); }				//��������ģ
} Vector;
double cross(const Point &a, const Point &b)					//���
{
	return  a.x*b.y - a.y*b.x;
}
double det(const Point &a, const Point &b)						//�����������ǲ�ˣ�det������ʽ����˼���������Ǹ�һ��
{
	return  a.x*b.y - a.y*b.x;
}
double dist(const Point &a, const Point &b)						//�����������
{
	return (a - b).norm();		//ֱ�ӷ���ģ�ͺ�
}
Vector rotate_Point(const Vector &p, double a)			//P������ʱ����תa(����)
{
	double tx = p.x, ty = p.y;							//t��ԭ��ֵ����������Ҹ�����p
	return Vector(tx*cos(a) - ty*sin(a), tx*sin(a) + ty*cos(a));		//x'=x*cos(a)-y*sin(a),y'=x*sin(a)+y*cos(a)
}

/*����Ϊֹ����������Ѿ��������*/

struct Line					//�����߶���
{
	Point a, b;
	Line(){}
	Line(Point x, Point y) :a(x), b(y){}
};
Line Point_make_Line(Point a, Point b)				//�������㣬����������ɵ��߶Ρ�
{
	return Line(a, b);
}
double dis_Point_segment(Point p, Point s, Point t)				//���p���߶�st�ľ���
{
	if (cmp((p - s)*(t - s)) < 0) return (p - s).norm();		//p�Ĵ�����st��s����ࣨsָ��p,t�����ֱ�߼нǴ���90�㣩
	if (cmp((p - t)*(s - t)) < 0) return (p - t).norm();		//p�Ĵ�����st��t����ࣨtָ��p,s�����ֱ�߼нǴ���90�㣩
	return fabs(cross(s - p, t - p) / dist(s, t));				//p�Ĵ�����st�ڣ��������������2�������߶γ�
}
double dis_Point_segment(Point p, Line l)				//���p���߶�l�ľ���,����һ��
{
	return dis_Point_segment(p, l.a, l.b);
}
Point Point_proj_Line(Point p, Point s, Point t)		//��p��ֱ��st�Ĵ��㣬��ȷ�������ǲ������߶��ϵ����·�Point_on_segment�ͺ�
{
	double r = (t - s)*(p - s) / ((t - s)* (t - s));			//sp��st�ϵ�ͶӰ���ȳ���st����
	return s + r*(t - s);								//s����ʼ�㣬����s->t�����ƶ���r�ĵ�
}
bool Point_on_segment(Point p, Point s, Point t)		//��p���߶�st�ϣ������˵㴦������true
{
	return (cmp(det(p - s, t - s)) == 0) && (cmp((p - s)*(p - t)) <= 0);	//���Ϊ0����p��s,t�м䣨�˵��ϣ�
}
bool parallel(Line a, Line b)							//��ֱ���Ƿ�ƽ�У�ֱ��ƽ�з���true
{
	return !cmp(det(a.a - a.b, b.a - b.b));				//a��������b��������˲�Ϊ0
}
bool Line_make_Point(Line a, Line b, Point &res)		//��ֱ���ཻ�Ľ��㣬�����ཻ����false����������res��
{
	if (parallel(a, b)) return false;
	double s1 = det(a.a - b.a, b.b - b.a);
	double s2 = det(a.b - b.a, b.b - b.a);
	res = (s1*a.b - s2*a.a) / (s1 - s2);				//���鿴����ѵĴ�׵� P257 ��ԭ��
	return true;
}
Line move_d(Line a, const double &len)					//��ֱ��a�ط��ߣ�a->b�����)����ƽ��len�ľ���õ��µ�ֱ��
{
	Vector d = a.b - a.a;
	d = d / d.norm();
	d = rotate_Point(d, Pi / 2.0);
	return Line(a.a + d*len, a.b + d*len);				//���ý����ˣ���㿴���ͺá�
}
/*���ˣ�������߶����Ѿ������ˣ��������Ƕ�����࣬��Ȼ�����̲���Ҫ����࣬������д��ȥ��*/
const int maxnn = 100;						//���������ж��ٸ����㣬�о�100�ھ������ǲ����õģ������ǿ��Ե�����
struct Polygon
{
	int n;								//�ö�����ж��ٸ�����
	Point a[maxnn];
	Polygon(){}
	double perimeter()					//�����ܳ�����
	{
		double sum = 0;					//�ܺ�
		a[n] = a[0];					//��0~n-1���㣬���԰�n����0����ȳɻ�
		for (int i = 0; i < n; i++) sum += (a[i + 1] - a[i]).norm();					//���߶ε�ģ����ȥ
		return sum;
	}
	double area()
	{
		double sum = 0;					//�ܺ�
		a[n] = a[0];					//ͬ��
		for (int i = 0; i < n; i++) sum += det(a[i + 1], a[i]);			//�ò�����������Ϊ������εĲ��ֲ��С��0��Ҳ����
		return sum / 2.0;				//��˽����������������������Ե��Ӻ󶼶�һ��ϵ��2��������
	}		//������Ϳ��Կ�����ѵĴ��P285�����������û��ѡ������ѵİ��ӡ��ܴ�̶��ϵ�ԭ��������ѵİ����Ҳ���Ϥ��

	int Point_in(Point t)						//�жϵ��Ƿ��ڶ�����ڣ����Ӷ�O(n)
	{											//����t����Ҫ�жϵĵ�t
		int num = 0, d1, d2, k;					//���0����t�ڶ������
		a[n] = a[0];							//���1����t�ڶ������
		for (int i = 0; i < n; i++)				//���2����t�ڶ���α���
		{
			if (Point_on_segment(t, a[i], a[i + 1])) return 2;		//�ڱ��Ϻܺÿ���
			k = cmp(det(a[i + 1] - a[i], t - a[i]));				//���k>0��ôt��a[i]->a[i+1]���
			d1 = cmp(a[i].y - t.y);
			d2 = cmp(a[i + 1].y - t.y);								//Ϊ�˺��Ե��������߼����d1,d2
			if (k > 0 && d1 <= 0 && d2 > 0) num++;
			if (k < 0 && d2 <= 0 && d1 > 0) num--;					//ÿ���߶���Ϊ���ұ�
		}
		return num != 0;
		/*�жϵ��ڶ�����ڣ��Ӹõ���һ��ˮƽ���ҵ����ߣ�ͬ�������������ཻ�����������ཻ����Ϊż����˵���������⣬�������ڡ�*/
		/*Ϊ�˱����غ��жϣ�������ϵ��߶���Ϊ���ֱ�*/
		/*ʵ���Ϻ������ֵ����⣬�ȷ�˵������һ����ͷ��״��������ʱ�����׳����⣬���Լ���/��ȥһ��΢Ԫ��>1e-8������������*/
	}
};
/*����εĺܶຯ���ҾͲ������ģ��������ˣ�ֻ�и��������࣬����Ϊֹ����β���Ҳ�������ˣ���������Բ*/

//�ܱ��ҵ��ǣ�Բ�ⶫ����Ҫ�Ĳ�����ȷ���࣬һ��Բ�ģ�һ���뾶/ֱ���� ���о�û�нṹ���ʾ��
//ԭ��Ʒ�е�dcmp��ͬ��cmp
//dot����dot
//abs����.norm,����ģ������˼
//Crosspet�����������ߵĽ���  ����ο�����135
//mysqrt�ǰ�0���µĶ�������0��Ū������
double mysqrt(double n)
{
	return sqrt(max(0.0, n));
}
int cirecle_cross_Line(Point a, Point b, Point o, double r, Point ret[])
{
	/*a,b�����߶�a->b��o��Բ�ģ�r�ǰ뾶��ret�Ǽ�������Ľ��㣬�����������У�����ֵ���ж��ٸ�����*/
	//|a+t(b-a)-o|=r Ȼ������ƽ���ͺã������һ��������ģ����ֱ��ƽ���������ǽ�������x����ƽ������y������ƽ��
	// sqr(x1 + t(x2-x1) - x0) + sqr( y1 + t(y2-y1) - y0) = r*r
	// sqr( (x1-x0)+t(x2-x1) ) == sqr(x1-x0) + 2*(x1-x0)*(x2-x1)*t + (x2-x1)*t*t
	//	�������x2-x1����dx��ͬ����y�Ϳ��Եõ�ϵ��������ϵ����������
	int num = 0;
	double x0 = o.x, y0 = o.y;
	double x1 = a.x, y1 = a.y;							//ֱ�߻��ɲ�����������aΪ���p=a+t*(b-a)��t�ǲ���
	double x2 = b.x, y2 = b.y;
	double dx = x2 - x1, dy = y2 - y1;					//Vector a->b���� (dx,dy)
	double A = dx*dx + dy*dy;							//a->b��ģ��ƽ��
	double B = 2 * dx * (x1 - x0) + 2 * dy * (y1 - y0);	//(o->a) * (a->b)*2
	double C = sqr(x1 - x0) + sqr(y1 - y0) - sqr(r);	//|o->a|ƽ��-r*r
	double delta = B * B - 4 * A * C;					//���η����Ƿ��и���A*t*t+B*t+C=0

	if (cmp(delta) >= 0)								//����ʵ������ȷ���м���
	{
		double t1 = (-B - mysqrt(delta)) / (2 * A);		//����t1��t2�ľ���ֵ�����ݹ�ʽ������
		double t2 = (-B + mysqrt(delta)) / (2 * A);
		if (cmp(t1 - 1) <= 0 && cmp(t1) >= 0)			//t1���ڣ�ͬʱt<=1��˵�����߶��ڣ���һ������
		{
			ret[num++] = Point(x1 + t1*dx, y1 + t1*dy);		//��������
		}
		if (cmp(t2 - 1) <= 0 && cmp(t2) >= 0)			//ͬ��t2���ڣ����߶���
		{
			ret[num++] = Point(x1 + t2*dx, y1 + t2*dy);		//��������
		}
		/*����ֻ��Ҫ��cmp(t1-1)<=0����cmp(t2-1)<=0ȥ���Ϳ����ж�ֱ����Բ�Ľ�������*/
	}
	return num;
}

struct Arc
{
	Point o;									//Բ��
	double r;									//�뾶
	double b, e;								//��ʼ���ȵ���ֹ���ȣ����մ���0�ļ������
	/*���ȵ��ص㣺�Ҷ�δ0����ʱ����ת��С��2��*/
	Arc() {}
	Arc(Point k, double f, double be, double ee) :o(k), r(f), b(be), e(ee){}	//���ع���

	double Point_in_Cirecle(Point k)			//����Բ�ڻ�����Բ��
	{
		if ((k - o).norm() <= r) return 1;		//��������ò�������
		return 0;
	}
	
	double angle(Point k)						//����o->k��x��Ļ���
	{
		Point b = k - o, a(1, 0);				//b��o->k,a��һ����׼��λ����
		double cosx = a*b / (a.norm()*b.norm());	//����ƫת����x��������һ��cos(x)
		double sinx = cross(a, b) / (a.norm()*b.norm());	//a��b�����b��a���������������Ǹ�
		if (sinx >= 0 && cosx >= 0)					//���Եĵ�һ����
			return acos(cosx);
		else if (sinx >= 0 && cosx < 0)				//�ڶ�����
			return acos(cosx);
		else if (sinx < 0 && cosx < 0)				//�������ޣ���ʱsinxС��0����ô�����æм���sinx�ķ����Ǻ����ľ���ֵ
			return Pi - asin(sinx);
		else										//�������ޣ�����ֵ����С��0��ֵ�õ�����
			return Pi * 2 + asin(sinx);
	}

	double dis(Point k)							//���һ���㵽һ��������̾���
	{
		double x = angle(k);
		if ((cmp(x - b) >= 0 && cmp(x - e) <= 0 && cmp(b - e) <= 0) || ((cmp(b - e) >= 0 && (cmp(x - b) >= 0 || cmp(x - e) <= 0))))//���߻ᾭ��Բ��
			//����ʼ�Ƕ�С����ֹ�Ƕ�ʱ��xӦ���������Ƕ��м�
			//����ʼ�Ƕȴ�����ֹ�Ƕȣ�˵��������x�������ᣬxӦ���������ᵽe������b��������֮��
		{
			return fabs(r - (o - k).norm());		//o->k��ģ�Ͱ뾶�Ĳ�ֵ
		}
		else										//������Բ������ȡһ
		{
			Point pb(cos(b)*r, sin(b)*r);			//��ʼ����Բ�����
			Point pe(cos(e)*r, sin(e)*r);			// ����ʱ���Բ���߶�
			return min((k - (o + pb)).norm(), (k - (o + pe)).norm());
		}
	}
};