#include<iostream>
#include<cmath>


using namespace std;


typedef long double LD;


template <typename T> class Punkt // punkt 2D
{
	private:
		T x;
		T y;
		
	public:
		Punkt(T __x = 0, T __y = 0) : x(__x), y(__y) {}
		
		~Punkt() {}
		
		Punkt(const Punkt <T> &other)
		{
			x = other.x;
			y = other.y;
		}
		
		bool operator == (Punkt <T> &other) const
		{
			return ( x == other.x  &&  y == other.y );
		}
		
		LD dist(Punkt <T> &other) const
		{
			LD dx = (x - other.x);
			LD dy = (y - other.y);
			
			return sqrt( dx*dx + dy*dy );
		}
		
		template <typename TT> friend LD dist(Punkt <TT> &A, Punkt <TT> &B);
		
		T get_x() const
		{
			return x;
		}
		
		T get_y() const
		{
			return y;
		}
		
		void set_x(T __x)
		{
			x = __x;
		}
		
		void set_y(T __y)
		{
			y = __y;
		}
		
		void set(T __x, T __y)
		{
			x = (__x);
			y = (__y);
		}
};


template <typename T> LD odl(Punkt <T> &A, Punkt <T> &B)
{
	LD dx = A.get_x() - B.get_x();
	LD dy = A.get_y() - B.get_y();
	
	return sqrt( dx*dx + dy*dy );
}


template <typename T> class Prostokat :public Punkt <T>
{
	private:
		int ID;
		Punkt <T> P;
		T dl;
		T wys;
		
	public:
		Prostokat()
		{
			ID = 0;
			
			P.set( 0, 0 );
			
			dl = 0;
			wys = 0;
		}
		
		Prostokat(int __ID, Punkt <T> &__P, T __dl, T __wys)
		{
			ID = ( __ID );
			P = ( __P );
			dl = ( __dl );
			wys = ( __wys );
		}
		
		~Prostokat() {}
		
		Prostokat(const Prostokat <T> &other) :Punkt<T> ()
		{
			ID = other.ID;
			P = other.P;
			dl = other.dl;
			wys = other.wys;
		}
		
		int get_ID() const
		{
			return ID;
		}
		
		void set_ID(int __ID)
		{
			ID = __ID;
		}
		
		Punkt <T> get_P() const
		{
			return P;
		}
		
		void set_P(Punkt <T> &__P)
		{
			P.set(  __P.get_x(), __P.get_y()  );
		}
		
		T get_dl() const
		{
			return dl;
		}
		
		void set_dl(T __dl)
		{
			dl = __dl;
		}
		
		T get_wys() const
		{
			return wys;
		}
		
		void set_wys(T __wys)
		{
			wys = __wys;
		}
		
		void set(int __ID, Punkt <T> &__P, T __dl, T __wys)
		{
			ID = ( __ID );
			P = ( __P );
			dl = ( __dl );
			wys = ( __wys );
		}
		
		T pole() const
		{
			return dl*wys;
		}
		
		T obw() const
		{
			return (T)(2) * ( dl + wys );
		}
		
		bool operator ==(Prostokat &other) const
		{
			return ( P == other.P  &&  dl == other.dl  &&  wys == other.wys );
		}
		
		T przeciecie(Prostokat &other) const
		{
			T x_lewy = max( P.get_x(), other.P.get_x() );
			T x_prawy = min( P.get_x() + dl, other.P.get_x() + other.dl );
			
			T dl = max( (T)(0), x_prawy - x_lewy );
			
			
			int y_dol = max( P.get_y(), other.P.get_y() );
			int y_gora = min( P.get_y() + wys, other.P.get_y() + other.wys );
			
			int wys = max( 0, y_gora - y_dol );
			
			// jesli czesc wspolna to prostokat o dadatnim polu, to lewy dolny wierzcholek ma wspolrzedne ( x_lewy, y_dol )
			
			return dl*wys;
		}
		
		LD przekatna() const
		{
			return sqrt( dl*dl + wys*wys );
		}
};



template <typename T> class Trojkat :public Punkt <T>
{
	private:
		Punkt <T> A;
		Punkt <T> B;
		Punkt <T> C;
		
	public:
		Trojkat()
		{
			A.set( 0, 0 );
			B.set( 0, 0 );
			C.set( 0, 0 );
		}
		
		Trojkat(Punkt <T> &__A, Punkt <T> &__B, Punkt <T> &__C)
		{
			A = __A;
			B = __B;
			C = __C;
		}
		
		~Trojkat() {}
		
		Trojkat(const Trojkat <T> &other) : Punkt <T> ()
		{
			A = other.A;
			B = other.B;
			C = other.C;
		}
		
		void set(Punkt <T> &__A, Punkt <T> &__B, Punkt <T> &__C)
		{
			A = __A;
			B = __B;
			C = __C;
		}
		
		LD obw() const
		{
			LD a = odl( A, B );
			LD b = odl( B, C );
			LD c = odl( C, A );
			
			return a+b+c;
		}
		
		LD pole() const
		{
			LD a = odl( A, B );
			LD b = odl( B, C );
			LD c = odl( C, A );
			
			LD p = (a+b+c) / 2. ;
			
			return sqrt(p) * sqrt(p-a) * sqrt(p-b) * sqrt(p-c);
		}
		
		LD hAB() const
		{
			LD P = (*this).pole();
			
			LD a = odl( A, B );
		
		
			return 2. * P / a;
		}
		
		LD hBC() const
		{
			LD P = (*this).pole();
			
			LD b = odl( B, C );
		
		
			return 2. * P / b;
		}
		
		LD hAC() const
		{
			LD P = (*this).pole();
			
			LD c = odl( A, C );
		
		
			return 2. * P / c;
		}
		
		template <typename TT> friend Punkt <LD> srodek_ciezkosci( Trojkat <T> &Triangle );
};


template <typename T> Punkt <LD> srodek_ciezkosci( Trojkat <T> &Triangle )
{
	Punkt <LD> S;
	
	LD xS = ( Triangle.A.get_x() + Triangle.B.get_x() + Triangle.C.get_x() ) / 3.;
	LD yS = ( Triangle.A.get_y() + Triangle.B.get_y() + Triangle.C.get_y() ) / 3.;
	
	S.set( xS, yS );
	
	
	return S;
}


template <typename T> class Wektor
{
	private:
		T x;
		T y;
	
	public:
		Wektor(T __x = 0, T __y = 0) : x(__x), y(__y) {}
		
		~Wektor() {}
		
		Wektor(const Wektor <T> &other)
		{
			x = other.x;
			y = other.y;
		}
		
		Punkt <T> get_x() const
		{
			return x;
		}
		
		Punkt <T> get_y() const
		{
			return y;
		}
		
		void set_x(T &__x)
		{
			x = __x;
		}
		
		void set_y(T &__y)
		{
			y = __y;
		}
		
		void set(T __x, T __y)
		{
			x = __x;
			y = __y;
		}
		
		LD dlugosc() const
		{
			return sqrt( x*x + y*y );
		}
		
		bool operator ==(Wektor <T> &other) const
		{
			return ( x == other.x  &&  y == other.y );
		}
		
		Wektor <T> operator +(Wektor <T> &other) const
		{
			Wektor <T> wyn;
			
			wyn.set( (*this).x + other.x, (*this).y + other.y );
			
			return wyn;
		}
		
		Wektor <T> operator -(Wektor <T> &other) const {
			Wektor <T> wyn;
			
			wyn.set( (*this).x - other.x, (*this).y - other.y );
			
			return wyn;
		}
		
		Wektor <T> operator +=(Wektor <T> &other)
		{
			Wektor <T> wyn;
			
			wyn.set( (*this).x + other.x, (*this).y + other.y );
			
			(*this) = wyn;
			
			return (*this);
		}
		
		Wektor <T> operator -=(Wektor <T> &other)
		{
			Wektor <T> wyn;
			
			wyn.set( (*this).x - other.x, (*this).y - other.y );
			
			(*this) = wyn;
			
			return (*this);
		}
		
		Wektor <T> operator *(T k) const
		{
			Wektor <T> wyn;
			
			wyn.set( k * (*this).x, k * (*this).y );
			
			return wyn;
		}
		
		Wektor <T> operator *=(T k)
		{
			Wektor <T> wyn;
			
			wyn.set( k * (*this).x, k * (*this).y );
			
			(*this) = wyn;
			
			return (*this);
		}
		
		
};


template <typename T> class Zesp
{
	private:
		T x;
		T y;
		
	public:
		Zesp(T  __x = 0, T  __y = 0) : x(__x), y(__y) {}
		
		~Zesp() {}
		
		Zesp(const Zesp <T> &other)
		{
			x = other.x;
			y = other.y;
		}
		
		T get_x() const
		{
			return x;
		}
		
		T get_y() const
		{
			return y;
		}
		
		void set_x(T __x)
		{
			x = __x;
		}
		
		void set_y(T __y)
		{
			y = __y;
		}
		
		void set(T __x, T __y)
		{
			x = __x;
			y = __y;
		}
		
		LD r() const
		{
			return sqrt( x*x + y*y );
		}
		
		bool operator ==(Zesp <T> &other) const
		{
			return ( x == other.x  &&  y == other.y );
		}
		
		Zesp <T> operator +(Zesp <T> &other) const
		{
			Zesp <T> wyn( (*this).x + other.x, (*this).y + other.y );
			
			return wyn;
		}
		
		Zesp <T> operator +=(Zesp <T> &other)
		{
			Zesp <T> wyn( (*this).x + other.x, (*this).y + other.y );
			
			(*this) = wyn;
			
			return (*this);
		}
		
		Zesp <T> operator -(Zesp <T> &other) const
		{
			Zesp <T> wyn( (*this).x - other.x, (*this).y - other.y );
			
			return wyn;
		}
		
		Zesp <T> operator -=(Zesp <T> &other)
		{
			Zesp <T> wyn( (*this).x - other.x, (*this).y - other.y );
			
			(*this) = wyn;
			
			return (*this);
		}
		
		Zesp <T> operator *(Zesp <T> &other) const
		{
			Zesp <T> wyn;
			
			wyn.set_x( x*other.x - y*other.y );
			wyn.set_y( y*other.x + x*other.y );
			
			return wyn;
		}
		
		Zesp <T> operator *=(Zesp <T> &other)
		{
			Zesp <T> wyn;
			
			wyn.set_x( x*other.x - y*other.y );
			wyn.set_y( y*other.x + x*other.y );
			
			
			(*this) = wyn;
			
			
			return (*this);
		}
		
		Zesp <T> sprzezenie() const
		{
			Zesp <T> wyn( x, -y );
			
			return wyn;
		}
};


int main()
{
	
	Punkt <int> P;
	
	P.set(9,7);
	
	cout << endl << P.get_x() << " " << P.get_y() << endl;
	
	
	Prostokat <int> R;
	
	R.set( 1, P, 8, 4 );
	
	cout << R.get_ID() << " " << R.get_dl() << " " << R.get_wys() << endl;
	
	return 0;
}
