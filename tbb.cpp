/*Two ways to compile
 * CC -dynamic  -std=c++0x -o fuck  tbb.cpp -tbb
 * or
 * module load tbb
 * export CRAYPE_LINK_TYPE=dynamic
 * CC tbb.cpp -tbb
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <vector>
#include <iostream>
#include <tbb/task_group.h>
#include <tbb/concurrent_unordered_map.h>

using namespace tbb;
using namespace std;

vector<int> vec;

void pushing(int arg){
    for(int i = 0; i <10; i++){
        vec.push_back(1);
    }
    cout<<vec.size()<<endl;

}


void cumpush(concurrent_unordered_multimap<int,int>& cum, int arg){
    for(int i = 0;i < 100; i++){
        cum.insert(make_pair(i%5,1));
    }
}

int Fib(int n) {
    if( n<2 ) {
            return n;
    } else {
            int x, y;
            tbb::task_group g;
            g.run([&]{x=Fib(n-1);}); // spawn a task
            g.run([&]{y=Fib(n-2);}); // spawn another task
            g.wait();                // wait for both tasks to complete
            cout<<x+y<<endl;
            return x+y;
            }
    }


int main(){
		concurrent_unordered_multimap<int,int> cum;
    pthread_t tid1, tid2;
   /* pthread_create(&tid1,NULL,pushing,NULL);
    pthread_create(&tid2,NULL,pushing,NULL);
    pthread_join(tid1,NULL);
    pthread_join(tid2,NULL);*/
    //Fib(3);
    task_group g;
    
    // g.run([&]{pushing((1));});
   // g.run([&]{pushing((1));});
    g.run([&]{
				cumpush(cum, 1);
		});
    
		g.run([&]{
					cumpush(cum, 1);
				});
 

    g.wait();
    cout<<cum.size()<<"count:"<<cum.count(0)<<endl;
    exit(0);

}
