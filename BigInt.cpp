#include <bits/stdc++.h>
using namespace std;
 

const int LengthOfGroup = 9; 
// group of digit working at a time


// Val_of_one_group=pow(10,LengthOfGroup)
const int Val_of_one_group = 1000000000;




#define ll long long int 

class bigint {
	//Storing data 
	// data stored in reverse order
public:
	vector<int> a;
	
	// Sign of a number 
	// num>=0  sign=   1
	// num<0   sign= (-1)
	int sign;
	

	// Length of number
	int Length(){
		if(a.size()==0) return 0;
		int len=(a.size()-1)*LengthOfGroup;
		int last=a[a.size()-1];
		while(last)
		{
			len++;
			last/=10;
		}
		return len;
	}


	//constructor
	bigint() :
		sign(1) {
	}
	
	//constructor
	bigint(long long v) {
		*this = v;
	}
	//constructor
	bigint(const string &s) {
		ReadInt(s);
	}

	//Assignment operation from bigInt
	void operator=(const bigint &v) {
		sign = v.sign;
		a = v.a;
	}
	

	//Assignment operation from long long
	void operator=(long long v) {
		sign = 1;
		a.clear();
		if (v < 0)
			sign = -1, v = -v;

		// long long can have more than Val_of_one_group value
		for (; v > 0; v = v / Val_of_one_group)
			a.push_back(v % Val_of_one_group);
	}
	

	// Readting and Converting string to parts of integer
	// stored in reverse order
	void ReadInt(const string &s) {
		sign = 1;
		a.clear();
		int pos = 0;
		if (pos < (int) s.size() && (s[pos] == '-' || s[pos] == '+')) {
			if (s[pos] == '-')
				sign = -1;
			pos=1;
		}
		
		for (int i = s.size() - 1; i >= pos; i -= LengthOfGroup) {
			int x = 0;
			for (int j = max(pos, i - LengthOfGroup + 1); j <= i; j++)
				x = x * 10 + s[j] - '0';
			a.push_back(x);
		}
		RemoveAllZeroFromBack();
	}
	
	// Overriding >> operator
	friend istream& operator>>(istream &stream, bigint &v) {
		string s;
		stream >> s;
		v.ReadInt(s);
		return stream;
	}
	
	// Overriding << operator
	// careful with print zero(0)
	friend ostream& operator<<(ostream &stream, const bigint &v) {
		if (v.sign == -1)
			stream << '-';
		stream << (v.a.empty() ? 0 : v.a.back());

		for (int i = (int) v.a.size() - 2; i >= 0; --i)
		{	
			if(v.a[i]!=0)
			stream <<  v.a[i];
			else    // if val=0 then print 0*9
			{
				for(int i=1;i<=LengthOfGroup;i++)
				stream<<"0";
			}
		}

		return stream;
	}

	// less than
	bool operator<(const bigint &v) const {
		if (sign != v.sign)
		{	// signs are diffrent then posative is bigger
			return sign < v.sign;
		}
		
		if (a.size() != v.a.size())
		{	// multiply with sign because of negative values
			return a.size() * sign < v.a.size() * v.sign;
		}

		for (int i = a.size() - 1; i >= 0; i--)
			if (a[i] != v.a[i])
			{	// multiply with sign because of negative values
				return a[i] * sign < v.a[i] * sign;
			}
		return false;
	}
	

	// reverse of less than
	bool operator>(const bigint &v) const {
		return v < *this;
	}


	bool operator<=(const bigint &v) const {
		return !(v < *this);
	}
	bool operator>=(const bigint &v) const {
		return !(*this < v);
	}
	bool operator==(const bigint &v) const {
		return !(*this < v) && !(v < *this);
	}
	bool operator!=(const bigint &v) const {
		return *this < v || v < *this;
	}

	// addition 
	bigint operator+(const bigint &v) const {
		if (sign == v.sign) {
			bigint res = v;
			for (int i = 0, carry = 0; i < (int) max(a.size(), v.a.size()) || carry; ++i) {
				if (i == (int) res.a.size())
				{	
					// adding one extra gap for operation in vector
					res.a.push_back(0);
				}
				res.a[i] += carry + (i < (int) a.size() ? a[i] : 0);
				carry = res.a[i] >= Val_of_one_group;
				if (carry)
				{
					res.a[i] -= Val_of_one_group;
				}
			}
			return res;
		}
		else
		{   // if sign are diff then addition is equal to subtration
			return *this - (-v);
		}
	}
	
	// subtration
	bigint operator-(const bigint &v) const {
		if (sign == v.sign) {
			if (abs() >= v.abs()) {
				bigint res = *this;
				for (int i = 0, carry = 0; i < (int) v.a.size() || carry; ++i) {
					res.a[i] -= carry + (i < (int) v.a.size() ? v.a[i] : 0);
					carry = res.a[i] < 0;
					if (carry)
						res.a[i] += Val_of_one_group;
				}
				res.RemoveAllZeroFromBack();
				return res;
			}
			return -(v - *this);
		}
		else
		{	// if sign are same then subtration is equal to addition
			return *this + (-v);
		}
	}



	// modulas by int val
	int operator%(int v) const {
			if (v < 0)
			{
				v = -v;
			}
		
		if(v==0){
		cout<<"mod by zero is not define"<<endl;
		assert(false);
		return 1;
		}

		int m = 0;
		for (int i = a.size() - 1; i >= 0; --i)
			m = (a[i] + m * (long long) Val_of_one_group) % v;
		return m * sign;
	}
 
	void operator+=(const bigint &v) {
		*this = *this + v;
	}
	void operator-=(const bigint &v) {
		*this = *this - v;
	}
	void operator*=(const bigint &v) {
		*this = *this * v;
	}
	void operator/=(const bigint &v) {
		*this = *this / v;
	}
	

	
 
	
 
	// chechking for zero value
	bool isZero() const {
		return a.empty() || (a.size() == 1 && !a[0]);
	}
	
	// making negative
	bigint operator-() const {
		bigint res = *this;
		res.sign = -sign;
		return res;
	}
	

	// absolute value
	// removing (-) sign
	bigint abs() const {
		bigint res = *this;
		res.sign =1;
		return res;
	}
	
	
	// Eucalideineuclidean algorithm
	friend bigint gcd(const bigint &a, const bigint &b) {
		return b.isZero() ? a : gcd(b, a % b);
	}
	friend bigint lcm(const bigint &a, const bigint &b) {
		return a / gcd(a, b) * b;
	}


	// fast power using binary exponantation
	bigint operator ^(const bigint &v){
		bigint ans=1,a=*this,b=v;
		while(!b.isZero()){
			if(b%2)
				ans*=a;
			a*=a,b/=2;
		}
		return ans;
	}


	// multiply with int
	void operator*=(int v) {
		if (v < 0)
			sign = -sign, v = -v;
		for (int i = 0, carry = 0; i < (int) a.size() || carry; ++i) {
			if (i == (int) a.size())
				a.push_back(0);
			long long cur = a[i] * (long long) v + carry;
			carry = (int) (cur / Val_of_one_group);
			a[i] = (int) (cur % Val_of_one_group);
		}
		RemoveAllZeroFromBack();
	}

	// multiply with int
	bigint operator*(int v) const {
		bigint res = *this;
		res *= v;
		return res;
	}


	// multiply with long long
	void operator*=(long long v) {
		if (v < 0)
			sign = -sign, v = -v;
		if(v > Val_of_one_group){
			*this = *this * (v / Val_of_one_group) * Val_of_one_group + *this * (v % Val_of_one_group);
			return ;
		}
		for (int i = 0, carry = 0; i < (int) a.size() || carry; ++i) {
			if (i == (int) a.size())
				a.push_back(0);
			long long cur = a[i] * (long long) v + carry;
			carry = (int) (cur / Val_of_one_group);
			a[i] = (int) (cur % Val_of_one_group);
		}
		RemoveAllZeroFromBack();
	}
	
	// multiply with long long
	bigint operator*(long long v) const {
		bigint res = *this;
		res *= v;
		return res;
	}


	// divion 
	bigint operator/(const bigint &v) const {
		return DivisonAndRem(*this, v).first;
	}
	//mod operration
	bigint operator%(const bigint &v) const {
		return DivisonAndRem(*this, v).second;
	}
	
	//divison with int  
	// beacause number is stored in reverse order
	void operator/=(int v) {
		if (v < 0)
			sign = -sign, v = -v;
		
		// i=a.Length()-1    => number stored in reverse order
		for (int i = (int) a.size() - 1, rem = 0; i >= 0; --i) {
			long long cur = a[i] + rem * (long long) Val_of_one_group;
			a[i] = (int) (cur / v);
			rem = (int) (cur % v);
		}
		RemoveAllZeroFromBack();
	}
	
	//divison with int
	bigint operator/(int v) const {
		bigint res = *this;
		res /= v;
		return res;
	}

	static vector<int> convert_Val_of_one_group(const vector<int> &a, int old_digits, int new_digits) {
		vector<long long> p(max(old_digits, new_digits) + 1);
		p[0] = 1;
		for (int i = 1; i < (int) p.size(); i++)
			p[i] = p[i - 1] * 10;
		vector<int> res;
		long long cur = 0;
		int cur_digits = 0;
		for (int i = 0; i < (int) a.size(); i++) {
			cur += a[i] * p[cur_digits];
			cur_digits += old_digits;
			while (cur_digits >= new_digits) {
				res.push_back(int(cur % p[new_digits]));
				cur /= p[new_digits];
				cur_digits -= new_digits;
			}
		}
		res.push_back((int) cur);
		while (!res.empty() && !res.back())
			res.pop_back();
		return res;
	}
 
	typedef vector<long long> vll;
	
	// Time complexity =  T(n)=3*T(n/2)+ O(n^1)
	// from master theoram Time complexity =O(n^(log 3 Val_of_one_group 2))  =O(n^(1.59))

	// c=1 a=3 b=2
	// log2(3)>c   ==> Time complexity =O(n^(log 3 Val_of_one_group 2))



	//  a= a1a2
	//  b= b1b2
    //  s1=a1b1 
	//  s2=a2b2
	//	s3=(a1+a2)*(b1+b2)
	//  s4=s3-s1-s2
    //  ans=s1*(10^(n))+s4*(10^(n/2))+s2
	//  n=length of a or b
	static vll MultipyWithKaratusba(const vll &a, const vll &b) {
		int n = a.size();
		vll res(n + n);
		
		if (n <= 16) // multiplication with normal method
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					res[i + j] += a[i] * b[j];
			return res;
		}
 
		int k = (n/2);
		vll a1(a.begin(), a.begin() + k);
		vll a2(a.begin() + k, a.end());
		vll b1(b.begin(), b.begin() + k);
		vll b2(b.begin() + k, b.end());
 
		vll a1b1 = MultipyWithKaratusba(a1, b1);
		vll a2b2 = MultipyWithKaratusba(a2, b2);
 
		for (int i = 0; i < k; i++)
			a2[i] += a1[i];
		for (int i = 0; i < k; i++)
			b2[i] += b1[i];
 
		vll r = MultipyWithKaratusba(a2, b2);
		for (int i = 0; i < (int) a1b1.size(); i++)
			r[i] -= a1b1[i];
		for (int i = 0; i < (int) a2b2.size(); i++)
			r[i] -= a2b2[i];
 
		for (int i = 0; i < (int) r.size(); i++)
			res[i + k] += r[i];
		for (int i = 0; i < (int) a1b1.size(); i++)
			res[i] += a1b1[i];
		for (int i = 0; i < (int) a2b2.size(); i++)
			res[i + n] += a2b2[i];
		return res;
	}
 
	bigint operator*(const bigint &v) const {
		
		vector<int> a6 = convert_Val_of_one_group(this->a, LengthOfGroup, 6);
		
		vector<int> b6 = convert_Val_of_one_group(v.a, LengthOfGroup, 6);
		
		vll a(a6.begin(), a6.end());
		vll b(b6.begin(), b6.end());
		while (a.size() < b.size())
			a.push_back(0);
		while (b.size() < a.size())
			b.push_back(0);

		// let a.size()=n
		//	n&(n-1) is the number removing last set bit from n
		//	(n&(n-1))==0  => number is power of 2 
		// and this is very useful in karatusba algo becase it work only for even numbers
		while (a.size() & (a.size() - 1))    // nearest power of 2 
			a.push_back(0), b.push_back(0);
		vll c = MultipyWithKaratusba(a, b);
		bigint res;
		res.sign = sign * v.sign;
		for (int i = 0, carry = 0; i < (int) c.size(); i++) {
			long long cur = c[i] + carry;
			res.a.push_back((int) (cur % 1000000));
			carry = (int) (cur / 1000000);
		}
		res.a = convert_Val_of_one_group(res.a, 6, LengthOfGroup);
		// removes Extra zeros from 
		res.RemoveAllZeroFromBack();
		return res;
	}

	void RemoveAllZeroFromBack() {
		while (!a.empty() && !a.back())
			a.pop_back();
		if (a.empty())
			sign = 1;
	}

	friend pair<bigint, bigint> DivisonAndRem(const bigint &a1, const bigint &b1) {
		
		
		
		bigint a = a1.abs();
		bigint b = b1.abs();
		bigint zero("0");
		if(b==zero) 
		{
		    cout<<"DIVISON BY ZERO ERROR\n";
		    assert(false);
		    
		}
		bigint q, r;
		q.a.resize(a.a.size());
	
	    // simple divison method 

		for (int i = a.a.size() - 1; i >= 0; i--) {
			r *= Val_of_one_group;
			r += a.a[i];
			int s1 = (r.a.size() <= b.a.size()) ? 0 : r.a[b.a.size()];
			int s2 = (r.a.size() <= b.a.size() - 1) ? 0 : r.a[b.a.size() - 1];

			int d = ((long long) Val_of_one_group * s1 + s2) / b.a.back();
			// d= approx q
			r -= b * d;
			while (r < 0)
			{
				r += b; 
				--d;
				cout<<"OK\n";
			}
			
			q.a[i] = d;
		}
 
		q.sign = a1.sign * b1.sign;
		r.sign = a1.sign;
		q.RemoveAllZeroFromBack();
		r.RemoveAllZeroFromBack();
		return make_pair(q, r);
	}
	
	
};


int main(){
	bigint a("125000");
	bigint b("0");

	
	cout<<"this is "<<(a/b)<<'\n';
}