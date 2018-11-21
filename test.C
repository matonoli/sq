using namespace std;

int wtf(int why) {
	return why;
}

int* func(int kof = 1) {

	cout << kof << endl;
	static int bla[3];
	bla[0] = kof*4;
	bla[1] = kof*5;
	bla[2] = kof*6;
	 //= {kof*4,kof*5,kof*6};
	return bla;
}

int test() {
	int* pbla = func(1);
	printf("pbla is %i , %i , %i \n", *(pbla),*(pbla+1),*(pbla+2));

	pbla = func(2);
	printf("pbla is %i , %i , %i \n", *(pbla),*(pbla+1),*(pbla+2));


	pbla = func(1);
	printf("pbla is %i , %i , %i \n", *(pbla),*(pbla+1),*(pbla+2));

}