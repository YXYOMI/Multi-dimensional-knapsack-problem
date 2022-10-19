//
//
//  Created by Xinyi OUYANG, 20030865
//  
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include<sys/time.h>

/* global parameters */
int RAND_SEED[] = {1,20,30,40,50,60,70,80,90,100,110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
int NUM_OF_RUNS = 1;
double MAX_TIME = 300;  //max amount of time permited (in sec)
int num_of_problems;

/* parameters for evlutionary algorithms*/
static int POP_SIZE = 50;   //please modify these parameters according to your problem
int MAX_NUM_OF_GEN = 100000; //max number of generations
float CROSSOVER_RATE = 0.9;
float MUTATION_RATE = 0.01;
float POP_RAND_RATE = 0.4;
int K= 3; // k-opt is used

/* declare parameters for variable neighbourhood search here*/
struct solution_struct best_sln;  //global best solution
struct solution_struct* population;  //global population
struct solution_struct* offspring;  //global offsprings of curent population

//return a random number between 0 and 1
float rand_01()
{
    float number;
    number = (float) rand();
    number = number/RAND_MAX;
    return number;
}

//return a random nunber ranging from min to max (inclusive)
int rand_int(int min, int max)
{
    int div = max-min+1;
    int val =rand() % div + min;
    return val;
}


struct item_struct{
    int dim; //no. of dimensions
    int* size; //volume of item in all dimensions
    int p;
    double ratio;
    int indx;
};

struct problem_struct{
    int n; //number of items
    int dim; //number of dimensions
    struct item_struct* items;
    int* capacities;  //knapsack capacities
};

void free_problem(struct problem_struct* prob)
{
    if(prob!=NULL)
    {
        if(prob->capacities !=NULL) free(prob->capacities);
        if(prob->items!=NULL)
        {
            for(int j=0; j<prob->n; j++)
            {
                if(prob->items[j].size != NULL)
                    free(prob->items[j].size);
            }
            free(prob->items);
        }
        free(prob);
    }
}

void init_problem(int n, int dim, struct problem_struct** my_prob)
{
    struct problem_struct* new_prob = malloc(sizeof(struct problem_struct));
    new_prob->n=n; new_prob->dim=dim;
    new_prob->items=malloc(sizeof(struct item_struct)*n);
    for(int j=0; j<n; j++)
        new_prob->items[j].size= malloc(sizeof(int)*dim);
    new_prob->capacities = malloc(sizeof(int)*dim);
    *my_prob = new_prob;
}


//example to create problem instances, actual date should come from file
struct problem_struct** load_problems(char* data_file)
{
    int i,j,k;
    //int num_of_probs;
    FILE* pfile = fopen(data_file, "r");
    if(pfile==NULL)
        {printf("Data file %s does not exist. Please check!\n", data_file); exit(2); }
    fscanf (pfile, "%d", &num_of_problems);
 
    struct problem_struct** my_problems = malloc(sizeof(struct problem_struct*)*num_of_problems);
    for(k=0; k<num_of_problems; k++)
    {
        int n, dim, obj_opt;
        fscanf (pfile, "%d", &n);
        fscanf (pfile, "%d", &dim); fscanf (pfile, "%d", &obj_opt);
        init_problem(n, dim, &my_problems[k]);  //allocate data memory
        for(j=0; j<n; j++)
        {
            my_problems[k]->items[j].dim=dim;
            my_problems[k]->items[j].indx=j;
            fscanf(pfile, "%d", &(my_problems[k]->items[j].p)); //read profit data
        }
        for(i=0; i<dim; i++)
        {
            for(j=0; j<n; j++)
            {
                fscanf(pfile, "%d", &(my_problems[k]->items[j].size[i])); //read size data
            }
        }
        for(i=0; i<dim; i++){
            fscanf(pfile, "%d", &(my_problems[k]->capacities[i]));
        }
    }
    fclose(pfile); //close file
    return my_problems;
}

struct solution_struct{
    struct problem_struct* prob; //maintain a shallow copy of the problem data
    float objective;
    int feasibility; //indicate the feasiblity of the solution
    int* x; //chromosome vector
    int* cap_left; //capacity left in all dimensions
    
};

void free_solution(struct solution_struct* sln)
{
    if(sln!=NULL)
    {
        free(sln->x);
        free(sln->cap_left);
        sln->objective=0;
        sln->prob=NULL;
        sln->feasibility=false;
    }
}


//copy a solution from another solution
bool copy_solution(struct solution_struct* dest_sln, struct solution_struct* source_sln)
{

    if(source_sln ==NULL) return false;
    if(dest_sln==NULL)
    {
        dest_sln = malloc(sizeof(struct solution_struct));
    }
    else{
        free(dest_sln->cap_left);
        free(dest_sln->x);
    }

    int n = source_sln->prob->n;
    int m =source_sln->prob->dim;
    dest_sln->x = malloc(sizeof(int)*n);
    dest_sln->cap_left=malloc(sizeof(int)*m);
    
    for(int i=0; i<m; i++)
        dest_sln->cap_left[i]= source_sln->cap_left[i];
    for(int j=0; j<n; j++)
        dest_sln->x[j] = source_sln->x[j];
    dest_sln->prob= source_sln->prob;
    dest_sln->feasibility=source_sln->feasibility;
    dest_sln->objective=source_sln->objective;
    return true;
}

void evaluate_solution(struct solution_struct* sln)
{
    //evaluate the feasibility and objective of the solution
    sln->objective =0; sln->feasibility = 1;
    struct item_struct* items_p = sln->prob->items;

    for(int i=0; i< items_p->dim; i++)
    {
        sln->cap_left[i]=sln->prob->capacities[i];
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->cap_left[i] -= items_p[j].size[i]*sln->x[j];
            if(sln->cap_left[i]<0) {
                sln->feasibility = -1*i; //exceeding capacity
                return;
            }
        }
    }

    if(sln->feasibility>0)
    {
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->objective += sln->x[j] * items_p[j].p;
        }
    }
}

//output a given solution to a file
void output_solution(struct solution_struct* sln, char* out_file)
{
    if(out_file !=NULL){
        FILE* pfile = fopen(out_file, "a"); //append solution data
        fprintf(pfile, "%i\n", (int)sln->objective);
        for(int i=0; i<sln->prob->n; i++)
        {
            fprintf(pfile, "%i ", sln->x[i]);
        }
        fprintf(pfile, "\n");
        /*for(int j=0; j<sln->prob->n; j++)
            fprintf(pfile, "%i ", sln->prob->items[j].p);
        fprintf(pfile, "\n");*/
        fclose(pfile);
    }
    else
        printf("sln.feas=%d, sln.obj=%f\n", sln->feasibility, sln->objective);
}


//check the  feasiblity and obj values of solutons from solution_file.
//return 0 is all correct or the index of the first infeasible problem [1, num_of_problems].
int check_solutions(struct problem_struct** my_problems, char* solution_file)
{
    FILE * pfile= fopen(solution_file, "r");
    if(pfile==NULL)
    {
        printf("Solution file %s does not exist. Please check!\n", solution_file);
        exit(2);
    }
    float val_obj;
    int val;
    fscanf (pfile, "%i", &val);
    if(val != num_of_problems)
    {
        printf("The stated number of solutions does not match the number of problems.\n");
        exit(3);
    }
    struct solution_struct temp_sln;
    int count=0, k=0;
    int n, dim;
    while(fscanf (pfile, "%f", &val_obj)!=EOF && k<num_of_problems)
    {
        //val_obj = val;
        n= my_problems[k]->n;  dim= my_problems[k]->dim;
        temp_sln.x = malloc(sizeof(int)*n);
        temp_sln.cap_left=malloc(sizeof(int)*dim);
        temp_sln.prob = my_problems[k];
        while(fscanf (pfile, "%i", &val)!=EOF)
        {
            if(val<0 || val>1) {fclose(pfile);  return k+1;} //illeagal values
            temp_sln.x[count] = val;
            count++;
            if(count==n)
            {
                evaluate_solution(&temp_sln);
                if(!temp_sln.feasibility || fabs(temp_sln.objective - val_obj)>0.01)
                {
                    fclose(pfile);
                    //printf("feasb=%i, obj= %f, val=%i\n",temp_sln.feasibility, temp_sln.objective, val_obj);
                    //output_solution(&temp_sln, "my_debug.txt");
                    return k+1;  //infeasible soltuion or wrong obj
                }
                else{
                    break;
                }
            }
        }
        count=0; k++;
        
        free(temp_sln.x); free(temp_sln.cap_left);
    }
    fclose(pfile);
    return 0;
}

//update global best solution from sln
void update_best_solution(struct solution_struct* sln)
{
    if(best_sln.objective < sln->objective)
    copy_solution(&best_sln, sln);
}

//check whether the two 1,0 items can swap
bool can_swap(struct solution_struct* sln, int item1, int item2)
{
    for(int d =0; d<sln->prob->dim; d++)
    {
        if(sln->cap_left[d]+sln->prob->items[item1].size[d] < sln->prob->items[item2].size[d])
            return false;
    }
    return true;
}


bool can_move(int nb_indx, int* move, struct solution_struct* curt_sln )
{
    bool ret=true;
    if(nb_indx==1)
    {
        int i = move[0];
        if(i<0) return false;
        for(int d=0; d<curt_sln->prob->dim; d++){
            if(curt_sln->cap_left[d] < curt_sln->prob->items[i].size[d])
                return false;
        }
    }
    else if(nb_indx==2){
        ret=can_swap(curt_sln, move[0], move[1]);
    }
    else if(nb_indx==3){//3-item swap
        int i= move[0], j= move[1], k= move[2];
        if(i<0 || j<0 || k<0) return false;
        if(curt_sln->x[j]>0) {//2-1 swap
            for(int d=0; d<curt_sln->prob->dim; d++){
                if(curt_sln->cap_left[d] + curt_sln->prob->items[i].size[d] +
                   curt_sln->prob->items[j].size[d] < curt_sln->prob->items[k].size[d])
                    return false;
            }
        }
        else {//1-2 swap
            for(int d=0; d<curt_sln->prob->dim; d++){
                if(curt_sln->cap_left[d] + curt_sln->prob->items[i].size[d] <
                   curt_sln->prob->items[j].size[d] + curt_sln->prob->items[k].size[d])
                    return false;
            }
        }
        
    }
    else ret=false;
    return ret;
}

bool apply_move(int nb_indx, int* move, struct solution_struct* sln )
{
    bool ret=true;
    if(nb_indx==1)
    {
        int i = move[0];
        if(i<0) return false;
        for(int d=0; d<sln->prob->dim; d++){
            sln->cap_left[d] -= sln->prob->items[i].size[d];
        }
        sln->objective += sln->prob->items[i].p;
        sln->x[i]=1;
        
        //printf("success\n");
    }
    else if(nb_indx==2){
        for(int d=0; d<sln->prob->dim; d++){
            sln->cap_left[d] = sln->cap_left[d] + sln->prob->items[move[0]].size[d]-
                sln->prob->items[move[1]].size[d];
        }
        sln->objective += sln->prob->items[move[1]].p-sln->prob->items[move[0]].p;
        sln->x[move[0]]=0; sln->x[move[1]]=1;
    }
    else if(nb_indx==3){//3-item swap
        int i= move[0], j= move[1], k= move[2];
        if(i<0 || j<0 || k<0) return false;
        if(sln->x[j]>0) {//2-1 swap
            for(int d=0; d<sln->prob->dim; d++){
                sln->cap_left[d] = sln->cap_left[d]+sln->prob->items[i].size[d] +
                    sln->prob->items[j].size[d] - sln->prob->items[k].size[d];
            }
            sln->objective += sln->prob->items[k].p-sln->prob->items[i].p-sln->prob->items[j].p;
            sln->x[i]=0; sln->x[j]=0; sln->x[k]=1;
        }
        else {//1-2 swap
            for(int d=0; d<sln->prob->dim; d++){
                sln->cap_left[d] = sln->cap_left[d]+sln->prob->items[i].size[d] -
                    sln->prob->items[j].size[d] - sln->prob->items[k].size[d];
            }
            sln->objective += sln->prob->items[j].p+sln->prob->items[k].p-sln->prob->items[i].p;
            sln->x[i]=0; sln->x[j]=1; sln->x[k]=1;
        }
        
    }
    else ret=false;
    return ret;
}

//nb_indx <=3
struct solution_struct* best_descent_vns(int nb_indx, struct solution_struct* curt_sln)
{
    struct solution_struct* best_neighb = malloc(sizeof(struct solution_struct));
    best_neighb->cap_left = malloc(sizeof(int)*curt_sln->prob->dim);
    best_neighb->x = malloc(sizeof(int)*curt_sln->prob->n);
    copy_solution(best_neighb, curt_sln);
    int n=curt_sln->prob->n;
    int curt_move[] ={-1,-1,-1}, best_move []={-1,-1,-1}, delta=0, best_delta=0;  //storing best neighbourhood moves
    int rand;
    switch (nb_indx)
    {
        case 1: //check whether any items can be inserted.
            for(int i=0; i<n; i++){
                if(curt_sln->x[i]>0) continue;
                curt_move[0]=i;
                if(can_move(nb_indx, &curt_move[0], best_neighb)){
                    delta = curt_sln->prob->items[i].p;
                    if(delta> best_delta) {
                        best_delta = delta; best_move[0] = i;
                    }
                }
            }
            if(best_delta>0) {    apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
        case 2:
            for(int i=0; i<n; i++){
                if(curt_sln->x[i]<=0) continue;
                for(int j=0; j<n; j++){
                    if(curt_sln->x[j]==0)
                    {
                        curt_move[0]= i; curt_move[1]= j; curt_move[2]=-1;
                        if(can_move(nb_indx, &curt_move[0], best_neighb)){
                            delta = curt_sln->prob->items[j].p -curt_sln->prob->items[i].p;
                            if(delta > best_delta){
                                best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=-1;
                            }
                        }
                    }
                }
            }
            if(best_delta>0) { apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
        case 3:
        //rand = rand_01();
        //if(rand < 0.5){
                //2-1 swap
            for(int i=0; i<n; i++){
                if(curt_sln->x[i]==0) continue;
                for(int j=0; j!=i&&j<n; j++){
                    if(curt_sln->x[j]==0) continue;
                    for(int k=0;k<n;k++){
                        if(curt_sln->x[k] == 0)
                        {
                            curt_move[0]=i; curt_move[1]=j; curt_move[2]=k;
                            if(can_move(nb_indx, &curt_move[0], best_neighb)){
                                delta = curt_sln->prob->items[k].p -curt_sln->prob->items[i].p-curt_sln->prob->items[j].p;
                                if(delta > best_delta){
                                    best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=k;
                                }
                            }
                        }
                    }
                }
            }
            //if(best_delta>0) { apply_move(nb_indx, &best_move[0], best_neighb);}
            //break;
        //}
            
        //if(rand >= 0.5){
                //1-2 swap
            for(int i=0; i<n; i++){
                if(curt_sln->x[i]==0) continue;
                for(int j=0; j<n; j++){
                    if(curt_sln->x[j]>0) continue;
                    for(int k=0;k!=j&&k<n;k++){
                        if(curt_sln->x[k] == 0)
                        {
                            curt_move[0]=i; curt_move[1]=j; curt_move[2]=k;
                            if(can_move(nb_indx, &curt_move[0], curt_sln)){
                                delta = curt_sln->prob->items[k].p +curt_sln->prob->items[j].p-curt_sln->prob->items[i].p;
                                if(delta > best_delta){
                                    best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=k;
                                }
                            }
                        }
                    }
                }
            }
            if(best_delta>0) { apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
        //}
        default:
            printf("Neighbourhood index is out of the bounds, nothing is done!\n");
    }
            
    return best_neighb;
}

//VNS
int varaible_neighbourhood_search()
{
    //do VNS in every individuals of offspring[].
    for(int p=0; p<POP_SIZE; p++){
        int nb_indx = 0;
        struct solution_struct* curt_sln = &offspring[p];
        update_best_solution(curt_sln);

        while(nb_indx<K){
            struct solution_struct* neighb_s=best_descent_vns(nb_indx+1, curt_sln); //best solution in neighbourhood nb_indx
            if(neighb_s->objective > curt_sln->objective){
                copy_solution(curt_sln, neighb_s);
                nb_indx=1;
            }
            else nb_indx++;
            free_solution(neighb_s);free(neighb_s);
        }
    }
    return 0;
}

//sort by ratio, in decreasing order
int cmpfunc1(const void* a, const void* b){
    const struct item_struct* item1 = a;
    const struct item_struct* item2 = b;
    if(item1->ratio>item2->ratio) return -1;
    if(item1->ratio<item2->ratio) return 1;
    return 0;
}

//sort by item index, in increasing order
int cmpfunc2 (const void * a, const void * b) {
    const struct item_struct* item1 = a;
    const struct item_struct* item2 = b;
    if(item1->indx>item2->indx) return 1;
    if(item1->indx<item2->indx) return -1;
    return 0;
}

//sort by ratio, in increasing order
int cmpfunc3(const void* a, const void* b){
    const struct item_struct* item1 = a;
    const struct item_struct* item2 = b;
    if(item1->ratio>item2->ratio) return 1;
    if(item1->ratio<item2->ratio) return -1;
    return 0;
}

//sort by objective, in decreasing order
int cmpfunc_sln (const void * a, const void * b){
    const struct solution_struct* sln1 = a;
    const struct solution_struct* sln2 = b;
    if(sln1->objective > sln2 ->objective) return -1;
    if(sln1->objective < sln2 ->objective) return 1;
    return 0;
}

//repair the solution by removing the least favorable items
void feasibility_repair(struct problem_struct* prob){

    //caculate the ratio of each items
    for(int i=0; i<prob->n;i++){
        double avg_size=0;
        struct item_struct* item_i = &prob->items[i];
        for(int d=0; d< prob->dim; d++){
            avg_size += (double)item_i->size[d]/prob->capacities[d];
        }
        item_i->ratio = item_i->p/avg_size;
    }
    qsort(prob->items, prob->n, sizeof(struct item_struct), cmpfunc3); //sort the item by ratio,in increasing order

    int ratio_index[prob->n]; // store item index when the items sorted by ratio
    for(int i=0; i<prob->n; i++){
        ratio_index[i] = prob->items[i].indx;
    }

    qsort(prob->items, prob->n, sizeof(struct item_struct), cmpfunc2); //resort the item original order by sorting by index 
    
    for(int k=0; k<POP_SIZE; k++){
        struct solution_struct* pop = &offspring[k];
        if(pop->feasibility != 1){
            for(int i=0; i<prob->n; i++){
                if(pop->x[ratio_index[i]] == 1){
                    pop->x[ratio_index[i]] = 0;
                    evaluate_solution(pop);
                    if(pop->feasibility != 1){
                        continue;
                    }else{
                        break;
                    }
                }else
                    continue;
            }
        }
    }
}

//a greedy heuristic solution based on profit/volume ratio
struct solution_struct* greedy_heuristic(struct problem_struct* prob)
{
    for(int i=0; i<prob->n;i++){
        double avg_size=0;
        struct item_struct* item_i = &prob->items[i];
        for(int d=0; d< prob->dim; d++){
            avg_size += (double)item_i->size[d]/prob->capacities[d];
        }
        item_i->ratio = item_i->p/avg_size;
    }
    qsort(prob->items, prob->n, sizeof(struct item_struct), cmpfunc1);
    
    struct solution_struct* init_sln = malloc(sizeof(struct solution_struct));
    init_sln->prob=prob;    init_sln->objective =0;
    init_sln->x = malloc(sizeof(int)*prob->n);
    init_sln->cap_left = malloc(sizeof(int)*prob->dim);
    int* cap = malloc(sizeof(int)*prob->dim);
    int i=0, d=0;
    for(d=0; d<prob->dim; d++) cap[d]=0; //aggregated volume
    for(i=0; i<prob->n; i++)
    {
        struct item_struct* item_i = &prob->items[i];
        //printf("item[%d].ratio = %.3f\t",item_i->indx,prob->items[i].ratio);
        for(d=0; d<prob->dim; d++){
            if(cap[d] + item_i->size[d] > prob->capacities[d])
                break; //infeasible to pack this item, try next
        }
        if(d>=prob->dim){
            init_sln->x[item_i->indx] = 1;
            init_sln->objective += item_i->p;
            for(d=0; d<prob->dim; d++){
                cap[d] += item_i->size[d];
            }
            //printf("packing item %d\n", item_i->indx);
        }
        else init_sln->x[item_i->indx] =0;
    }
    for(d=0; d<prob->dim; d++){
        init_sln->cap_left[d] = prob->capacities[d]- cap[d];
    }
    free(cap);
    //restore item original order by sorting by index.
    qsort(prob->items, prob->n, sizeof(struct item_struct), cmpfunc2);
    
    evaluate_solution(init_sln);
    //output_solution(init_sln, "greedy_sln.txt");
    //printf("Init_sln obj=\t%f\tfeasiblity = %d.\n", init_sln->objective, init_sln->feasibility);
    return init_sln;
}

//tournament select
void tournament(){

    int tournament = 3;
    int rand[tournament];
    int temp;
    
    for(int k=0; k<POP_SIZE; k++){
        for(int i=0; i<tournament; i++)
            rand[i] = rand_int(0, POP_SIZE-1);
        for(int i=1; i<tournament-1; i++){
            if(population[rand[0]].objective < population[rand[i]].objective)
                rand[0] = rand[i];  //store the index of population with largest objective in rand[0]
        }
        copy_solution(&offspring[k], &population[rand[0]]); //the best individual copy to the mating pool
    }
}

//uniform crossover
void ux(struct solution_struct* parent1, struct solution_struct* parent2){
    
    int chorom_len = parent1->prob->n;
    float rand;
    int temp;

    for(int i=0; i<chorom_len; i++){
        rand = rand_01();
        if(rand >= 0.5){
            temp = parent1->x[i];
            parent1->x[i] = parent2->x[i];
            parent2->x[i] = temp;
        }
    }
    evaluate_solution(parent1);
    evaluate_solution(parent2);
}

void crossover(){

    float rand;
    int select = 0;
    int parent1, parent2;

    for(int i=0; i<POP_SIZE; i++){
        rand = rand_01();
        if(rand < CROSSOVER_RATE){
            select++;
            if((select%2) == 0){
                parent2 = i;  //find parent2 and do the uniform crossover
                ux(&offspring[parent1], &offspring[parent2]);
            }else if(i!= POP_SIZE-1)
                parent1 = i; //find parent1
        }
    }

}

void mutate(int chorom_len){
    
    for(int i=0; i<POP_SIZE; i++){
        for(int j=0; j<chorom_len; j++){
            float rand = rand_01();
            if(rand < MUTATION_RATE){
                offspring[i].x[j] = !(offspring[i].x[j]);
            }
        }
        evaluate_solution(&offspring[i]);
    }
}

//initial population by greedy and random generate, the size of random individuals is dependent on POP_RAND_RATE
void generate(struct problem_struct* prob){
    
    int pop_rand_size = (int)POP_SIZE * POP_RAND_RATE;
    int pop_neighb_size = POP_SIZE - pop_rand_size;

    struct solution_struct* curt_sln = greedy_heuristic(prob);
    update_best_solution(curt_sln);

    int i, j;
    int count = 0;
    
    //find the neighbours of current solution (from greedy search, 1-1 swap)
        for(i=0; i<prob->n; i++){
            if(count >= pop_neighb_size)
                break;
            if(count == 0){
                copy_solution(&population[0], curt_sln);
                count++;
            }else{
            if(curt_sln->x[i]>0){
                for(j=0; j<prob->n; j++){           
                    if(i!=j && curt_sln->x[j]==0 && can_swap(curt_sln,i, j))
                    {
                        copy_solution(&population[count], curt_sln);
                        for(int d=0; d<prob->dim; d++){
                            population[count].cap_left[d] = population[count].cap_left[d]+ prob->items[i].size[d]-prob->items[j].size[d];
                        }
                        int delta =curt_sln->prob->items[j].p -curt_sln->prob->items[i].p;
                        population[count].objective = curt_sln->objective + delta;
                        population[count].x[i] = 0; //swap
                        population[count].x[j] = 1;
                        count++;
                        if(count >= pop_neighb_size)
                            break;
                    }
                }
            }
            }
        }

    int can_change, rand, iter;
    
    //initial the rand individuals (all bits of chromosome are set to 0)
    for(i=count; i<POP_SIZE; i++){
        for(j=0; j<prob->n; j++){
            population[i].x[j] = 0;
        }
    }
    
    
    while(count < POP_SIZE){
        population[count].prob = prob;
        //initial the cap_left
        for(i=0; i<prob->dim; i++){
            population[count].cap_left[i] = prob->capacities[i];
        }
        
        //do the random generate
        for(iter=0; iter<1000; iter++){
            rand = rand_int(0, prob->n-1);
            if(population[count].x[rand] != 1){
                for(i=0; i<prob->dim; i++){
                    can_change = 0;
                    population[count].cap_left[i] -= prob->items[rand].size[i];
                    if(population[count].cap_left[i] < 0)
                        break;
                    can_change = 1;
                }
                if(can_change == 1)
                    population[count].x[rand] = 1;
            }else
                continue;
        }
        count++;
    }
        
    for(i=0; i<POP_SIZE; i++){
        evaluate_solution(&population[i]);
    }

    free_solution(curt_sln);
    free(curt_sln);

}

//find the best solution in the offspring
void find_best_sln(){
    for(int i=0; i<POP_SIZE; i++){
        update_best_solution(&population[i]);
    }
}

//combine current population and its offspring, and generate next population
void replace(struct problem_struct* prob){

    int i, j;
    struct solution_struct* par_and_son = malloc(sizeof(struct solution_struct) * (POP_SIZE*2));

    for(int k=0; k<POP_SIZE*2; k++)
    {
        par_and_son[k].cap_left= malloc(sizeof(int)*prob->dim);
        par_and_son[k].x = malloc(sizeof(int)*prob->n);
        par_and_son[k].objective = -10000;
    }

    for(i=0; i<POP_SIZE; i++){
        copy_solution(&par_and_son[i], &population[i]);
    }
        
    for(i=0; i<POP_SIZE; i++)
        copy_solution(&par_and_son[i+POP_SIZE], &offspring[i]);

    qsort(par_and_son, POP_SIZE*2, sizeof(struct solution_struct), cmpfunc_sln); //sort the parents and offspring by objective, in decreasing order

    for(i=0; i<POP_SIZE; i++)
        copy_solution(&population[i], &par_and_son[i]);
    
    for(int k=0; k<POP_SIZE*2; k++)
        free_solution(&par_and_son[k]);
    free(par_and_son);

}

int MemeticAlgorithm(struct problem_struct* prob){

    int gen = 0;
    clock_t time_start, time_fin;
    time_start = clock();
    double time_spent = 0;

    generate(prob);  //initial population
    while(gen < MAX_NUM_OF_GEN && time_spent < MAX_TIME){
        gen++;
        //printf("gegneration: %d, ", gen);
        tournament();  //do tournament select
        crossover();   //do uniform crossover
        mutate(prob->n);  //mutate
        feasibility_repair(prob);    //repair the offspring
        varaible_neighbourhood_search();  //do the local search (VNS)

        replace(prob);   //generate the next population
        find_best_sln(); //update best solution in the population
        
        time_fin=clock();
        time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;

        //printf("time: %f, ", time_spent);
        //printf("best: %f\n", best_sln.objective);
    }
    
    return 0;
}

int main(int argc, const char * argv[]) {
    printf("Starting the run!\n");
    char data_file[50]={"somefile"}, out_file[50]={}, solution_file[50]={};  //max 50 problem instances per run
    if(argc<3)
    {
        printf("Insufficient arguments. Please use the following options:\n   -s data_file (compulsory)\n   -o out_file (default my_solutions.txt)\n   -c solution_file_to_check\n   -t max_time (in sec)\n");
        return 1;
    }
    else if(argc>9)
    {
        printf("Too many arguments.\n");
        return 2;
    }
    else
    {
        for(int i=1; i<argc; i=i+2)
        {
            if(strcmp(argv[i],"-s")==0)
                strcpy(data_file, argv[i+1]);
            else if(strcmp(argv[i],"-o")==0)
                strcpy(out_file, argv[i+1]);
            else if(strcmp(argv[i],"-c")==0)
                strcpy(solution_file, argv[i+1]);
            else if(strcmp(argv[i],"-t")==0)
                MAX_TIME = atoi(argv[i+1]);
        }
        //printf("data_file= %s, output_file= %s, sln_file=%s, max_time=%d", data_file, out_file, solution_file, MAX_TIME);
    }
    struct problem_struct** my_problems = load_problems(data_file);

    if(strlen(solution_file)<=0)
    {
        if(strcmp(out_file,"")==0) strcpy(out_file, "my_solutions.txt"); //default output
        FILE* pfile = fopen(out_file, "w"); //open a new file
        fprintf(pfile, "%d\n", num_of_problems); fclose(pfile);
        
        for(int k=0; k<num_of_problems; k++)
        {
            population = malloc(sizeof(struct solution_struct) * POP_SIZE);
            offspring = malloc(sizeof(struct solution_struct) * POP_SIZE);
            for(int m=0; m<POP_SIZE; m++)
            {
                population[m].cap_left= malloc(sizeof(int)*my_problems[k]->dim);
                population[m].x = malloc(sizeof(int)*my_problems[k]->n);
                offspring[m].cap_left= malloc(sizeof(int)*my_problems[k]->dim);
                offspring[m].x = malloc(sizeof(int)*my_problems[k]->n);
            }
            printf("Problem[%d]\n", k);
            best_sln.objective=0; best_sln.feasibility=0;
            for(int run=0; run<NUM_OF_RUNS; run++)
            {
                srand(RAND_SEED[2]);
                MemeticAlgorithm(my_problems[k]); // call MA method
            }
            printf("best_sln: %f\n\n", best_sln.objective);
            output_solution(&best_sln,out_file);
            for(int k=0; k<POP_SIZE; k++)
                free_solution(&population[k]);
            free(population);

            for(int k=0; k<POP_SIZE; k++)
                free_solution(&offspring[k]);
            free(offspring);
        }
    }
    for(int k=0; k<num_of_problems; k++)
    {
       free_problem(my_problems[k]); //free problem data memory
    }
    free(my_problems); //free problems array

    

    if(best_sln.x!=NULL && best_sln.cap_left!=NULL){ free(best_sln.cap_left); free(best_sln.x);} //free global
    
    return 0;
}
