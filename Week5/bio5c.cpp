#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
using namespace std;


string input[21][21];

int table[128][128];

string s1,s2;


typedef pair<int,int> P1;
typedef pair<P1,P1> P2;

const int MAXN = 10010;
int dp[MAXN];
int dp2[MAXN];
int sumdp[MAXN];
const int indel = 5;


struct Solver{
	const string &s1;
	const string &s2;
	int n,m,mid;

	int f[2][MAXN];
	int b[2][MAXN];
	int fcur,bcur;

	Solver(const string &s1,const string &s2):s1(s1),s2(s2){
		n = s1.length();
		m = s2.length();
		mid = m/2;
		fcur = 0;
		bcur = 0;
	}
	void update_dp_forward(int j,const int *from,int *to){
		for(int i=0;i<n;i++) {
			to[i]=-10000000;
			to[i] = max(to[i],from[i]-indel);
			if(i>0 && j>0 ) to[i] = max(to[i],from[i-1]+ table[s1[i-1]][s2[j-1]]);
			if(i>0) to[i] = max(to[i],to[i-1]-indel);
		}
	}
	void update_dp_backward(int j,const int *from,int *to){
		for(int i=n-1;i>=0;i--) {
			to[i]=-10000000;
			to[i] = max(to[i],from[i]-indel);
			if(i+1<n && j<m ) to[i] = max(to[i],from[i+1]+ table[s1[i]][s2[j]]);
			if(i+1<n) to[i] = max(to[i],to[i+1]-indel);
		}
	}


	P1 middle_node(	){

		for(int i=0;i<n;i++) f[0][i]=b[0][n-i-1] = -i*indel;

		for(int i=0;i<=mid;i++) {
			update_dp_forward(i,f[fcur],f[fcur^1]);
			fcur^=1;
		}

		for(int i=m-1;i>mid;i--){
			update_dp_backward(i,b[bcur],b[bcur^1]);
			bcur^=1;
		}
		int tm = -1000000;
		int p=0;
		for(int i=0;i<m;i++) {
			if(f[fcur][i]+b[bcur][i]>tm) {
				tm = f[fcur][i]+b[bcur][i];
				p=i;
			}

		}
		return P1(p,mid);
	}



	P1 getnext(P1 node){
		int i = node.first;
		int j = node.second;
		if(table[s1[i]][s2[j]]> - indel){
			return P1(i+1,j+1);
		}else return P1(i+1,j);
	}
};



int match[MAXN];
void solve(string s1,string s2,int off1 = 0,int off2=0){

	if(s1.length()==0 || s2.length() ==0 ) return;

	Solver solver(s1,s2);
	P1 ans = solver.middle_node();
	match[ans.first+off1]=ans.second+off2;

	string sub1 = s1.substr(0,ans.first);
	string sub2 = s2.substr(0,ans.second);

	string sub3 = s1.substr(ans.first+1,string::npos);
	string sub4 = s2.substr(ans.second+1,string::npos);



	cout<<sub1<<":"<<sub2<<endl;
	cout<<sub3<<":"<<sub4<<endl;
	solve(sub1,sub2,off1,off2);
	solve(sub3,sub4,off1+ans.first+1,off2+ans.second+1);
}

int main() {
	setbuf(stdout,NULL);
	memset(table,0xf0,sizeof(table));
	for(int i=0;i<21;i++){
		for(int j=0;j<21;j++) cin>>input[i][j];
	}

	for(int i=1;i<=20;i++)
		for(int j=1;j<=20;j++){
			table[input[i][0][0]][input[0][j][0]] = atoi(input[i][j].c_str());
		}

	cin>>s1>>s2;


	cout<<"S1="<<s1<<endl;
	cout<<"S2="<<s2<<endl;
	Solver solver(s1,s2);

	solve(s1,s2);


	for(int i=0;i<s1.length();i++) printf("%d %d\n",i,match[i]);
	return 0;
}
