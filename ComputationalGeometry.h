double max(double a, double b)						//重载一次max
{
	if (a > b) return a;
	return b;
}
double min(double a, double b)
{
	if (a > b) return  b;
	return a;
}
const double Pi = acos(-1.0);						//Pi是π
const double eps = 1e-8;							//接近0的一个误差允许范围
int cmp(double x)									//计算几何误差修正函数
{
	if (fabs(x) < eps) return 0;
	if (x>0) return 1;
	return -1;
}
inline double sqr(double x) { return x*x; }			//平方
typedef struct  Point			//点坐标
{
	double x, y;
	Point(){}
	Point(double a, double b) :x(a), y(b){}
	void input(){ scanf("%lf%lf", &x, &y); }	//没什么用的读入
	/*对向量/点的加减，向量常数乘除，向量点乘，向量/点的比较进行了重载*/
	friend Point operator + (const Point &a, const Point &b){ return Point(a.x + b.x, a.y + b.y); }
	friend Point operator - (const Point &a, const Point &b){ return Point(a.x - b.x, a.y - b.y); }
	friend bool operator == (const Point &a, const Point &b){ return (!cmp(a.x - b.x)) && (!cmp(a.y - b.y)); }
	friend Point operator * (const Point &a, const double &b){ return Point(a.x * b, a.y * b); }
	friend Point operator * (const double &a, const Point &b){ return Point(a * b.x, a * b.y); }
	friend double operator * (const Point &a, const Point &b){ return (a.x * b.x + a.y * b.y); }
	friend Point operator / (const Point &a, const double &b){ return Point(a.x / b, a.y / b); }
	double norm(){ return sqrt(sqr(x) + sqr(y)); }				//这个是求个模
} Vector;
double cross(const Point &a, const Point &b)					//叉乘
{
	return  a.x*b.y - a.y*b.x;
}
double det(const Point &a, const Point &b)						//这两个函数是叉乘，det是行列式的意思，和上面那个一样
{
	return  a.x*b.y - a.y*b.x;
}
double dist(const Point &a, const Point &b)						//计算两点距离
{
	return (a - b).norm();		//直接返回模就好
}
Vector rotate_Point(const Vector &p, double a)			//P向量逆时针旋转a(弧度)
{
	double tx = p.x, ty = p.y;							//t是原有值，不能随便乱改引用p
	return Vector(tx*cos(a) - ty*sin(a), tx*sin(a) + ty*cos(a));		//x'=x*cos(a)-y*sin(a),y'=x*sin(a)+y*cos(a)
}

/*至此为止，点和向量已经定义完成*/

struct Line					//定义线段类
{
	Point a, b;
	Line(){}
	Line(Point x, Point y) :a(x), b(y){}
};
Line Point_make_Line(Point a, Point b)				//给两个点，返还它们组成的线段。
{
	return Line(a, b);
}
double dis_Point_segment(Point p, Point s, Point t)				//求点p到线段st的距离
{
	if (cmp((p - s)*(t - s)) < 0) return (p - s).norm();		//p的垂足在st的s的外侧（s指向p,t两点的直线夹角大于90°）
	if (cmp((p - t)*(s - t)) < 0) return (p - t).norm();		//p的垂足在st的t的外侧（t指向p,s两点的直线夹角大于90°）
	return fabs(cross(s - p, t - p) / dist(s, t));				//p的垂足在st内，用三角形面积的2倍除以线段长
}
double dis_Point_segment(Point p, Line l)				//求点p到线段l的距离,重载一次
{
	return dis_Point_segment(p, l.a, l.b);
}
Point Point_proj_Line(Point p, Point s, Point t)		//求p到直线st的垂足，想确定垂足是不是在线段上调用下方Point_on_segment就好
{
	double r = (t - s)*(p - s) / ((t - s)* (t - s));			//sp在st上的投影长度除以st长度
	return s + r*(t - s);								//s是起始点，沿着s->t方向移动了r的点
}
bool Point_on_segment(Point p, Point s, Point t)		//点p在线段st上（包括端点处）返回true
{
	return (cmp(det(p - s, t - s)) == 0) && (cmp((p - s)*(p - t)) <= 0);	//叉乘为0，且p在s,t中间（端点上）
}
bool parallel(Line a, Line b)							//两直线是否平行，直线平行返还true
{
	return !cmp(det(a.a - a.b, b.a - b.b));				//a的向量与b的向量叉乘不为0
}
bool Line_make_Point(Line a, Line b, Point &res)		//两直线相交的交点，若不相交返还false，交点存放在res中
{
	if (parallel(a, b)) return false;
	double s1 = det(a.a - b.a, b.b - b.a);
	double s2 = det(a.b - b.a, b.b - b.a);
	res = (s1*a.b - s2*a.a) / (s1 - s2);				//详情看刘汝佳的大白的 P257 有原理。
	return true;
}
Line move_d(Line a, const double &len)					//将直线a沿法线（a->b的左侧)方向平移len的距离得到新的直线
{
	Vector d = a.b - a.a;
	d = d / d.norm();
	d = rotate_Point(d, Pi / 2.0);
	return Line(a.a + d*len, a.b + d*len);				//懒得解释了，随便看看就好。
}
/*到此，红书的线段类已经结束了，接下来是多边形类，虽然本工程不需要这个类，但是先写上去吧*/
const int maxnn = 100;						//多边形最多有多少个顶点，感觉100在竞赛中是不够用的，所以是可以调整的
struct Polygon
{
	int n;								//该多边形有多少个顶点
	Point a[maxnn];
	Polygon(){}
	double perimeter()					//计算周长函数
	{
		double sum = 0;					//总和
		a[n] = a[0];					//从0~n-1个点，所以把n点与0点相等成环
		for (int i = 0; i < n; i++) sum += (a[i + 1] - a[i]).norm();					//把线段的模加上去
		return sum;
	}
	double area()
	{
		double sum = 0;					//总和
		a[n] = a[0];					//同上
		for (int i = 0; i < n; i++) sum += det(a[i + 1], a[i]);			//用叉乘求面积，因为凹多边形的部分叉乘小于0，也能用
		return sum / 2.0;				//叉乘结果是三角形面积两倍，线性叠加后都多一个系数2，除掉。
	}		//具体解释可以看刘汝佳的大白P285，但这个代码没有选用刘汝佳的板子。很大程度上的原因是刘汝佳的板子我不熟悉。

	int Point_in(Point t)						//判断点是否在多边形内，复杂度O(n)
	{											//输入t，需要判断的点t
		int num = 0, d1, d2, k;					//输出0代表t在多边形外
		a[n] = a[0];							//输出1代表t在多边形内
		for (int i = 0; i < n; i++)				//输出2代表t在多边形边上
		{
			if (Point_on_segment(t, a[i], a[i + 1])) return 2;		//在边上很好看懂
			k = cmp(det(a[i + 1] - a[i], t - a[i]));				//如果k>0那么t在a[i]->a[i+1]左侧
			d1 = cmp(a[i].y - t.y);
			d2 = cmp(a[i + 1].y - t.y);								//为了忽略掉左侧的射线加入的d1,d2
			if (k > 0 && d1 <= 0 && d2 > 0) num++;
			if (k < 0 && d2 <= 0 && d1 > 0) num--;					//每条线段设为左开右闭
		}
		return num != 0;
		/*判断点在多边形内：从该点做一条水平向右的射线，同级射线与多边形相交的情况，如果相交次数为偶数，说明点在形外，奇数在内。*/
		/*为了便于重合判断，多边形上的线段设为左开又闭*/
		/*实际上好像会出现点问题，比方说出现了一个箭头形状，过顶点时很容易出问题，可以加上/减去一个微元（>1e-8）解决这个问题*/
	}
};
/*多边形的很多函数我就不往这个模板里加入了，只有个基本的类，至此为止多边形部分也差不多结束了，接下来是圆*/

//很悲惨的是，圆这东西需要的参数的确不多，一个圆心，一个半径/直径， 书中就没有结构体表示。
//原作品中的dcmp等同于cmp
//dot就是dot
//abs就是.norm,就是模长的意思
//Crosspet好像是两条线的交点  详情参考红书135
//mysqrt是把0以下的东西返还0，弄掉虚数
double mysqrt(double n)
{
	return sqrt(max(0.0, n));
}
int cirecle_cross_Line(Point a, Point b, Point o, double r, Point ret[])
{
	/*a,b代表线段a->b，o是圆心，r是半径，ret是计算出来的交点，保存在数组中，返回值是有多少个交点*/
	//|a+t(b-a)-o|=r 然后两边平方就好，左边是一个向量的模不能直接平方向量，是将向量的x坐标平方加上y向量的平方
	// sqr(x1 + t(x2-x1) - x0) + sqr( y1 + t(y2-y1) - y0) = r*r
	// sqr( (x1-x0)+t(x2-x1) ) == sqr(x1-x0) + 2*(x1-x0)*(x2-x1)*t + (x2-x1)*t*t
	//	把上面的x2-x1带成dx，同处理y就可以得到系数，具体系数计算如下
	int num = 0;
	double x0 = o.x, y0 = o.y;
	double x1 = a.x, y1 = a.y;							//直线化成参数方程依照a为起点p=a+t*(b-a)，t是参数
	double x2 = b.x, y2 = b.y;
	double dx = x2 - x1, dy = y2 - y1;					//Vector a->b就是 (dx,dy)
	double A = dx*dx + dy*dy;							//a->b的模的平方
	double B = 2 * dx * (x1 - x0) + 2 * dy * (y1 - y0);	//(o->a) * (a->b)*2
	double C = sqr(x1 - x0) + sqr(y1 - y0) - sqr(r);	//|o->a|平方-r*r
	double delta = B * B - 4 * A * C;					//二次方程是否有根，A*t*t+B*t+C=0

	if (cmp(delta) >= 0)								//存在实根，不确定有几个
	{
		double t1 = (-B - mysqrt(delta)) / (2 * A);		//计算t1和t2的具体值，根据公式法运算
		double t2 = (-B + mysqrt(delta)) / (2 * A);
		if (cmp(t1 - 1) <= 0 && cmp(t1) >= 0)			//t1存在，同时t<=1，说明在线段内，有一个交点
		{
			ret[num++] = Point(x1 + t1*dx, y1 + t1*dy);		//存入数组
		}
		if (cmp(t2 - 1) <= 0 && cmp(t2) >= 0)			//同理t2存在，在线段内
		{
			ret[num++] = Point(x1 + t2*dx, y1 + t2*dy);		//存入数组
		}
		/*我们只需要把cmp(t1-1)<=0或者cmp(t2-1)<=0去掉就可以判断直线与圆的交点坐标*/
	}
	return num;
}

struct Arc
{
	Point o;									//圆心
	double r;									//半径
	double b, e;								//起始弧度到终止弧度，按照大于0的计算好了
	/*弧度的特点：右端未0，逆时针旋转，小于2π*/
	Arc() {}
	Arc(Point k, double f, double be, double ee) :o(k), r(f), b(be), e(ee){}	//重载构造

	double Point_in_Cirecle(Point k)			//点在圆内或者在圆上
	{
		if ((k - o).norm() <= r) return 1;		//好像这次用不到来着
		return 0;
	}
	
	double angle(Point k)						//计算o->k与x轴的弧度
	{
		Point b = k - o, a(1, 0);				//b是o->k,a是一个基准单位向量
		double cosx = a*b / (a.norm()*b.norm());	//假设偏转角是x，我们算一个cos(x)
		double sinx = cross(a, b) / (a.norm()*b.norm());	//a×b，如果b在a左侧就是正，否则是负
		if (sinx >= 0 && cosx >= 0)					//明显的第一象限
			return acos(cosx);
		else if (sinx >= 0 && cosx < 0)				//第二象限
			return acos(cosx);
		else if (sinx < 0 && cosx < 0)				//第三象限，此时sinx小于0，那么我们用π加上sinx的反三角函数的绝对值
			return Pi - asin(sinx);
		else										//第四象限，满的值加上小于0的值得到弧度
			return Pi * 2 + asin(sinx);
	}

	double dis(Point k)							//求出一个点到一个弧的最短距离
	{
		double x = angle(k);
		if ((cmp(x - b) >= 0 && cmp(x - e) <= 0 && cmp(b - e) <= 0) || ((cmp(b - e) >= 0 && (cmp(x - b) >= 0 || cmp(x - e) <= 0))))//连线会经过圆弧
			//当起始角度小于终止角度时，x应该在两个角度中间
			//当起始角度大于终止角度，说明经过了x轴正半轴，x应该在正半轴到e或者是b到正半轴之间
		{
			return fabs(r - (o - k).norm());		//o->k的模和半径的差值
		}
		else										//不经过圆弧，二取一
		{
			Point pb(cos(b)*r, sin(b)*r);			//开始处的圆弧起点
			Point pe(cos(e)*r, sin(e)*r);			// 结束时候的圆弧线段
			return min((k - (o + pb)).norm(), (k - (o + pe)).norm());
		}
	}
};