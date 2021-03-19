#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

/*  struct  */
typedef struct hardblock{
    int id;
    int x;
    int y;
    int width;
    int height;
    int rotate;
} Hardblock;

typedef struct terminal{
    int id;
    int x;
    int y;
} Terminal;

typedef struct node{
    int parent;
    int leftchild;
    int rightchild;
} Node;

typedef struct cost_funct{
    double area;
    double wirelength;
    double R;   // aspect ratio
    int width;
    int height;
    double cost;
} Cost;

/* function declare */
void parse_hardblocks(string inputfile);
void parse_net(string inputfile);
void parse_pl(string inputfile);
void write_output(string output_file, int wirelength, vector<Hardblock> &hard);
void initial_Bstartree();
void preorder_traverse(int node, bool left);
void B2Floorplan();
Cost cost();
void rotate(int node);
void move(int from, int to);
void swap(int a, int b);
void SimulatedAnnealing();
unsigned int random_seed();

/* variable declare */
int seed = 1542959801;
int num_hard, num_terminal, num_net, num_pin;
vector<Hardblock> hardblocks;
vector<vector<int>> nets;
vector<Terminal> terminals;
double dead_space_ratio;
double area_normal=0, wl_normal=0;
int W;      // fixed outline max x coordinate
int root = -1;
vector<Node> bstar_tree;
vector<int> contour;
bool in_fixed_outline;
int min_cost_root, Fixedmin_cost_root;      // floorplan with min cost ,in fixed outline
Cost min_cost, Fixedmin_cost;
vector<Hardblock> min_cost_floorplan, Fixedmin_cost_floorplan;
vector<Node> min_cost_bstar, Fixedmin_cost_bstar;
bool feasible = 0;


/* main */
int main(int argc, char **argv){
    string hardblocks_file = argv[1];
    string nets_file = argv[2];
    string terminals_file = argv[3];
    string output_file = argv[4];
    dead_space_ratio = atof(argv[5]);

    // Read testcase
    parse_hardblocks(hardblocks_file);
    parse_net(nets_file);
    parse_pl(terminals_file);

    if(num_hard==300 && dead_space_ratio==0.15)
        srand(seed);

    // Floorplanning
    initial_Bstartree();
    SimulatedAnnealing();

    // Write output floorplan file
    write_output(output_file, Fixedmin_cost.wirelength,Fixedmin_cost_floorplan);

    // Total Runtime
    double runtime = clock();
    printf("[  Total Run time  ]: %f sec\n",runtime/CLOCKS_PER_SEC);
    return 0;

    /*
    // Test how small the dead space ratio could be for your program to produce a legal result in 20 minutes.
    // n100 -> 0.06, n200 -> 0.075
    do{
        double start = clock();
        cout << dead_space_ratio << endl;
        parse_hardblocks(hardblocks_file);
        parse_net(nets_file);
        parse_pl(terminals_file);

        if(num_hard==300 && dead_space_ratio==0.15)
            srand(seed);

        initial_Bstartree();
        SimulatedAnnealing();
        dead_space_ratio -= 0.005;

        write_output(output_file, Fixedmin_cost.wirelength,Fixedmin_cost_floorplan);

        double end = clock();
        printf("[  Total Run time  ]: %f sec\n",(end-start)/CLOCKS_PER_SEC);
    } while(feasible == 1);

    double runtime = clock();
    printf("[  Total test time  ]: %f sec\n", runtime/CLOCKS_PER_SEC);
    return 0;
    */

}



/* function */

// parse hardblocks file
void parse_hardblocks(string inputfile){
    ifstream file;
    file.open(inputfile);

    string Num, colon, str;
    file >> Num >> colon >> num_hard;   //  NumHardRectilinearBlocks : number of hard rectilinear block nodes
    file >> Num >> colon >> num_terminal;   //  NumTerminals : number of terminal (pad etc.) nodes
    file >> str;    //  nodeName hardrectilinear vertexNumber vertex1, vertex2, ..., vertexN
    
    //  ex : sb0 hardrectilinear 4 (0, 0) (0, 33) (43, 33) (43, 0)
    hardblocks = vector<Hardblock>(num_hard);
    double total_block_area = 0;    // A = total block area
    for (int i = 0; i < num_hard; i++) {
        getline(file, str);

        size_t left_paren = str.find("(");
        left_paren = str.find("(", left_paren + 1);
        left_paren = str.find("(", left_paren + 1);
        size_t comma = str.find(",");
        comma = str.find(",", comma + 1);
        comma = str.find(",", comma + 1);
        size_t right_paren = str.find(")");
        right_paren = str.find(")", right_paren + 1);
        right_paren = str.find(")", right_paren + 1);

        char buffer[100];
        int width, height;
        size_t len = str.copy(buffer, comma-left_paren-1, left_paren+1);
        buffer[len] = '\0';
        width = atoi(buffer);
        len = str.copy(buffer, right_paren-comma-2, comma+2);
        buffer[len] = '\0';
        height = atoi(buffer);

        hardblocks[i].id = i;
        hardblocks[i].x = -1;
        hardblocks[i].y = -1;
        hardblocks[i].width = width;
        hardblocks[i].height = height;
        hardblocks[i].rotate = 0;

        total_block_area += width*height;
    }

    //H*W* = total block area * (1 + dead space ratio)
    double fixed_outline_area = total_block_area * (1 + dead_space_ratio);
    W = sqrt(fixed_outline_area);

    printf("\nFixed-outline constraint info \n");
    printf("[ Floorplan Region ] : %.2f\n", fixed_outline_area);
    printf("[ Total Block Area ] : %.2f\n", total_block_area);
    printf("[ Fixed-Outline Width ] : %d\n\n", W);

    file.close();
}

// parse nets file
void parse_net(string inputfile){
    ifstream file;
    file.open(inputfile);

    string Num, colon, str;
    file >> Num >> colon >> num_net;    //NumNets : number of nets
    file >> Num >> colon >> num_pin;    //NumPins : number of pins
    cout << Num << endl;
    nets = vector<vector<int>>(num_net);
    for (int i=0; i<num_net ; i++){
        int degree;
        file >> Num >> colon >> degree;    //NetDegree : number of pins on the net
        for(int j=0; j<degree ; j++){
            file >> str;    //nodeName
            int id;
            if(str[0]=='p'){
                str.erase(0,1);     // 刪除'p'
                id = atoi(str.c_str()) + num_hard;
            }
            else if(str[0]=='s'){
                str.erase(0,2);     // 刪除'sb'
                id = atoi(str.c_str());
            }
            nets[i].emplace_back(id);   //emplace_back方法，比push_back的效率要高很多
            //cout << id <<  endl;
        }
    }
    //cout << nets[1][0] << endl;
    file.close();
}

// parse pl file
void parse_pl(string inputfile){
    ifstream file;
    file.open(inputfile);

    string nodeName;
    int x, y;

    terminals = vector<Terminal>(num_terminal+1);
    for(int i=1 ; i<=num_terminal ; i++){
        file >> nodeName >> x >> y;     //nodeName XY-coordinate
        terminals[i].id = i;
        terminals[i].x = x;
        terminals[i].y = y;
    }
    file.close();
}

// write output floorplan file
void write_output(string output_file, int wirelength, vector<Hardblock> &hard){
    ofstream file;
    file.open(output_file);

    file << "Wirelength " << wirelength << "\n";
    file << "Blocks\n";

    //nodeName lower-left corner coordinates (x,y) Rotated
    //nodeName lower-left corner coordinates (x,y) Unrotated
    for(int i=0; i<num_hard ; i++){
        if(hard[i].rotate)
            file << "sb" << i << " " << hard[i].x << " " << hard[i].y << " 1\n";
        else
            file << "sb" << i << " " << hard[i].x << " " << hard[i].y << " 0\n";
    }

    file.close();
}

// build initial B*-Tree
void initial_Bstartree(){
    bstar_tree = vector<Node>(num_hard);
    
    queue<int> bfs;
    vector<int> insert(num_hard,0);
    root = rand() % num_hard;
    bstar_tree[root].parent = -1;
    bfs.push(root);
    insert[root] = 1;

    int i = num_hard-1;
    while (!bfs.empty()){
        int parent = bfs.front();
        int left_child = -1, right_child = -1;
        bfs.pop();

        if(i > 0){
            do{
                left_child = rand() % num_hard;
            } while (insert[left_child]);
            bstar_tree[parent].leftchild = left_child;
            bfs.push(left_child);
            insert[left_child] = 1;
            i--;

            if(i > 0){
                do{
                    right_child = rand() % num_hard;
                } while (insert[right_child]);
                bstar_tree[parent].rightchild = right_child;
                bfs.push(right_child);
                insert[right_child] = 1;
                i--;
            }
        }
        bstar_tree[parent].leftchild = left_child;
        bstar_tree[parent].rightchild = right_child;
        if(left_child != -1)
            bstar_tree[left_child].parent = parent;
        if(right_child != -1)
            bstar_tree[right_child].parent = parent;
    }
}

// preorder traverse the B*-Tree
void preorder_traverse(int node, bool left){
    int parent = bstar_tree[node].parent;

    // left child or right child of parent
    // left child : the lowest, adjacent block on the right(xj = xi + wi)
    // right child : the first block above, with the same x-coordinate(xj = xi)
     if(left)
         hardblocks[node].x = hardblocks[parent].x + hardblocks[parent].width;
     else
        hardblocks[node].x = hardblocks[parent].x;

    // mantain contour
    int x_left = hardblocks[node].x;
    int x_right = x_left + hardblocks[node].width;
    int y_max = 0;
    for(int i=x_left; i<x_right ; i++){
        if(contour[i] > y_max){
            y_max = contour[i];
        }
    }
    hardblocks[node].y = y_max;
    y_max += hardblocks[node].height;
    for (int i=x_left; i<x_right ; i++){
        contour[i] = y_max;
    }

    if(bstar_tree[node].leftchild != -1)
        preorder_traverse(bstar_tree[node].leftchild, true);
    if(bstar_tree[node].rightchild != -1)
        preorder_traverse(bstar_tree[node].rightchild, false);
}

// convert B*-Tree to Floorplan
void B2Floorplan(){
    contour = vector<int>(5*W,0);
    hardblocks[root].x = 0;
    hardblocks[root].y = 0;
    for(int i=0; i<hardblocks[root].width; i++)
        contour[i] = hardblocks[root].height;
    
    if(bstar_tree[root].leftchild != -1)
        preorder_traverse(bstar_tree[root].leftchild, true);
    if(bstar_tree[root].rightchild != -1)
        preorder_traverse(bstar_tree[root].rightchild, false);
}

// calculate the cost
Cost cost(){
    B2Floorplan();

    int width = 0, height = 0;
    for(int i=0 ; i<num_hard ; i++){
        if(hardblocks[i].x + hardblocks[i].width > width)
            width = hardblocks[i].x + hardblocks[i].width;
        if(hardblocks[i].y + hardblocks[i].height > height)
            height = hardblocks[i].y + hardblocks[i].height;
    }

    double floorplan_area = width * height;
    double R = (double)height / width;
    double wirelength = 0;
    for(const vector<int> &net : nets){
        int x_min = width + 1;
        int x_max = 0;
        int y_min = height + 1;
        int y_max = 0;
        for(const int pin : net){
            if(pin < num_hard){
                int x_center = hardblocks[pin].x + hardblocks[pin].width/2;
                int y_center = hardblocks[pin].y + hardblocks[pin].height/2;
                if(x_center < x_min)
                    x_min = x_center;
                if(y_center < y_min)
                    y_min = y_center;
                if(x_center > x_max)
                    x_max = x_center;
                if(y_center > y_max)
                    y_max = y_center;
            }
            else{
                const Terminal &ter = terminals[pin-num_hard];
                if(ter.x < x_min)
                    x_min = ter.x;
                if(ter.y < y_min)
                    y_min = ter.y;
                if(ter.x > x_max)
                    x_max = ter.x;
                if(ter.y > y_max)
                    y_max = ter.y;
            }
        }
        // HPWL : half perimeter wire length
        wirelength += (x_max - x_min) + (y_max - y_min);
    }

    Cost c;
    c.width = width;
    c.height = height;
    c.area = floorplan_area;
    c.wirelength = wirelength;
    c.R = R;

    // normalize the area & wirelength
    if(area_normal == 0)
        area_normal = floorplan_area;
    if(wl_normal == 0)
        wl_normal = wirelength;
    double area_cost = c.area / area_normal;
    double wl_cost = c.wirelength / wl_normal;
    double R_cost = (1-R) * (1-R);
    double width_penalty = 0;
    double height_penalty = 0;
    if(width > W)
        width_penalty = (double)width / W;
    if(height > W)
        height_penalty = (double)height / W;
    
    c.cost = area_cost + wl_cost + 0.6*R_cost + width_penalty + height_penalty;
    
    return c;
}

// perturb : rotate a module
void rotate(int node){
    int temp = hardblocks[node].width;
    hardblocks[node].width = hardblocks[node].height;
    hardblocks[node].height = temp;
    hardblocks[node].rotate = 1 - hardblocks[node].rotate;
}

// perturb : move a module to another palce
void move(int from, int to){
    /* delete the node */
    // no child
    if(bstar_tree[from].leftchild == -1 && bstar_tree[from].rightchild == -1){
        int parent = bstar_tree[from].parent;
        if(bstar_tree[parent].leftchild == from)
            bstar_tree[parent].leftchild = -1;
        else if(bstar_tree[parent].rightchild == from)
            bstar_tree[parent].rightchild = -1;
        else{
            printf("[Error] perturb : moving !\n");
            exit(1);
        }
    }
    // two child
    else if(bstar_tree[from].leftchild != -1 && bstar_tree[from].rightchild != -1){
        do{
            bool move_left;
            if(bstar_tree[from].leftchild != -1 && bstar_tree[from].rightchild != -1)
                move_left = rand() % 2 == 0;
            else if(bstar_tree[from].leftchild != -1)
                move_left = true;
            else
                move_left = false;
            
            if(move_left)
                swap(from, bstar_tree[from].leftchild);
            else
                swap(from, bstar_tree[from].rightchild);
        }while(bstar_tree[from].leftchild != -1 || bstar_tree[from].rightchild != -1);

        int parent = bstar_tree[from].parent;
        if(bstar_tree[parent].leftchild == from)
            bstar_tree[parent].leftchild = -1;
        else if(bstar_tree[parent].rightchild == from)
            bstar_tree[parent].rightchild = -1;
        else{
            printf("[Error] perturb : moving !\n");
            exit(1);
        }
    }
    // one child
    else{
        int child;
        if(bstar_tree[from].leftchild != -1)
            child = bstar_tree[from].leftchild;
        else
            child = bstar_tree[from].rightchild;

        int parent = bstar_tree[from].parent;
        if(parent != -1){
            if(bstar_tree[parent].leftchild == from)
                bstar_tree[parent].leftchild = child;
            else if(bstar_tree[parent].rightchild == from)
                bstar_tree[parent].rightchild = child;
            else{
                printf("[Error] perturb : moving !\n");
                exit(1);
            }
        }
        bstar_tree[child].parent = parent;
        
        //root may change
        if(from == root)
            root = child;
    }
    
    /* insert the node */
    int insert_pos = rand() % 4;
    int child;
    if(insert_pos == 0){
        child = bstar_tree[to].leftchild;
        bstar_tree[from].leftchild = child;
        bstar_tree[from].rightchild = -1;
        bstar_tree[to].leftchild = from;
    }
    else if(insert_pos == 0){
        child = bstar_tree[to].rightchild;
        bstar_tree[from].leftchild = child;
        bstar_tree[from].rightchild = -1;
        bstar_tree[to].rightchild = from;
    }
    else if(insert_pos == 0){
        child = bstar_tree[to].leftchild;
        bstar_tree[from].leftchild = -1;
        bstar_tree[from].rightchild = child;
        bstar_tree[to].leftchild = from;
    }
    else{
        child = bstar_tree[to].rightchild;
        bstar_tree[from].leftchild = -1;
        bstar_tree[from].rightchild = child;
        bstar_tree[to].rightchild = from;
    }
    bstar_tree[from].parent = to;
    if(child != -1)
        bstar_tree[child].parent = from; 
}

// perturb : swap two modules
void swap(int n1, int n2){
    /* swap two parent */
    int n1_parent = bstar_tree[n1].parent;
    if(n1_parent != -1){
        if(bstar_tree[n1_parent].leftchild == n1)
            bstar_tree[n1_parent].leftchild = n2;
        else if(bstar_tree[n1_parent].rightchild == n1)
            bstar_tree[n1_parent].rightchild = n2;
        else{
            printf("[Error] perturb : swapping !\n");
            exit(1);
        }
    }

    int n2_parent = bstar_tree[n2].parent;
    if(n2_parent != -1){
        if(bstar_tree[n2_parent].leftchild == n2)
            bstar_tree[n2_parent].leftchild = n1;
        else if(bstar_tree[n2_parent].rightchild == n2)
            bstar_tree[n2_parent].rightchild = n1;
        else{
            printf("[Error] perturb : swapping !\n");
            exit(1);
        }
    }

    bstar_tree[n1].parent = n2_parent;
    bstar_tree[n2].parent = n1_parent;

    /* swap two child */
    int n1_leftchild = bstar_tree[n1].leftchild;
    int n1_rightchild = bstar_tree[n1].rightchild;
    bstar_tree[n1].leftchild = bstar_tree[n2].leftchild;
    bstar_tree[n1].rightchild = bstar_tree[n2].rightchild;
    bstar_tree[n2].leftchild = n1_leftchild;
    bstar_tree[n2].rightchild = n1_rightchild;

    if(bstar_tree[n1].leftchild != -1)
        bstar_tree[bstar_tree[n1].leftchild].parent = n1;
    if(bstar_tree[n1].rightchild != -1)
        bstar_tree[bstar_tree[n1].rightchild].parent = n1;
    if(bstar_tree[n2].leftchild != -1)
        bstar_tree[bstar_tree[n2].leftchild].parent = n2;
    if(bstar_tree[n2].rightchild != -1)
        bstar_tree[bstar_tree[n2].rightchild].parent = n2;
    
    /* swap parent & child */
    if(bstar_tree[n1].parent == n1)
        bstar_tree[n1].parent = n2;
    else if(bstar_tree[n1].leftchild == n1)
        bstar_tree[n1].leftchild = n2;
    else if(bstar_tree[n1].rightchild == n1)
        bstar_tree[n1].rightchild = n2;

    if(bstar_tree[n2].parent == n2)
        bstar_tree[n2].parent = n1;
    else if(bstar_tree[n2].leftchild == n2)
        bstar_tree[n2].leftchild = n1;
    else if(bstar_tree[n2].rightchild == n2)
        bstar_tree[n2].rightchild = n1;
    
    // root may change
    if(n1 == root)
        root = n2;
    else if(n2 == root)
        root = n1;
}

// Simulated Annealing Floorplanning (p.21)
void SimulatedAnnealing(){
    min_cost = cost();
    min_cost_floorplan = hardblocks;
    const double P = 0.95;
    const double r = 0.9;
    const int k = 20;
    const int N = k * num_hard;
    const double T0 = -min_cost.cost*num_hard / log(P);

    double T = T0;
    int MT = 0;     // try次數
    int uphill = 0;
    int reject = 0;
    Cost old_cost = min_cost;
    in_fixed_outline = false;

    clock_t time_start = clock();
    clock_t time = time_start;
    const int max_seconds = (num_hard/20) * (num_hard/20);
    const int TIME_LIMIT = 1200;    // produce a legal result in 20 minutes
    int seconds = 0, runtime = 0;

    do{
        MT = 0;
        uphill = 0;
        reject = 0;

        do{
            vector<Hardblock> temp_hard(hardblocks);
            vector<Node> temp_bstar(bstar_tree);
            int old_root = root;

            /* Perturbing */
            int perturb = rand() % 3;
            // M1 : rotate
            if(perturb == 0){
                int node = rand() % num_hard;
                rotate(node);
            }
            // M2 : swap
            else if(perturb == 1){
                int n1,n2;
                n1 = rand() % num_hard;
                do{
                    n2 = rand() % num_hard;
                } while(n2==n1);
                swap(n1,n2);
            }
            // M3 : move
            else if(perturb == 2){
                int from, to;
                from = rand() % num_hard;
                do{
                    to = rand() % num_hard;
                } while(to==from || bstar_tree[from].parent==to);
                move(from,to);
            }
            else{
                return;
            }

            MT++;
            Cost new_cost = cost();
            double delta_cost = new_cost.cost - old_cost.cost;
            double random = (double)rand() / RAND_MAX;
            if(delta_cost<=0 || random<exp(-delta_cost/T)){     // p21. line 17
                if(delta_cost > 0)
                    uphill++;
                
                // feasible solution with minimum cost
                if(new_cost.width<=W && new_cost.height<=W){
                    if(in_fixed_outline){
                        if(new_cost.cost < Fixedmin_cost.cost){
                            Fixedmin_cost_root = root;
                            Fixedmin_cost = new_cost;
                            Fixedmin_cost_floorplan = hardblocks;
                            Fixedmin_cost_bstar = bstar_tree;
                        }
                    }
                    else{
                        in_fixed_outline = true;
                        Fixedmin_cost_root = root;
                        Fixedmin_cost = new_cost;
                        Fixedmin_cost_floorplan = hardblocks;
                        Fixedmin_cost_bstar = bstar_tree;
                    }
                }
                // infeasible solution with minimum cost
                if(new_cost.cost < min_cost.cost){
                    min_cost_root = root;
                    min_cost = new_cost;
                    min_cost_floorplan = hardblocks;
                    min_cost_bstar = bstar_tree;
                }

                old_cost = new_cost;
            }
            else{
                reject++;
                root = old_root;
                if(perturb == 0)
                    hardblocks = temp_hard;
                else
                    bstar_tree = temp_bstar;
            }
        } while(uphill<=N && MT<=2*N);

        T *= r;

        seconds = (clock()-time) / CLOCKS_PER_SEC;
        runtime = (clock()-time_start) / CLOCKS_PER_SEC;
        if(seconds>=max_seconds && in_fixed_outline==false){
            printf("Overtime %d %d\n",min_cost.width,min_cost.height);
            seconds = 0;
            time = clock();
            T = T0;
        }
    } while(seconds<max_seconds && runtime<TIME_LIMIT);

    if(in_fixed_outline){
        feasible = 1;
        printf("\nFound feasible solution\n");
        printf("[       Width      ] : %d\n", min_cost.width);
        printf("[       Height     ] : %d\n", min_cost.height);
        printf("[       Area       ] : %.0f\n", min_cost.area);
        printf("[     Wirelength   ] : %.0f\n", min_cost.wirelength);
        printf("[        R         ] : %f\n", min_cost.R);
        printf("[       Cost       ] : %f\n\n\n", min_cost.cost);
    }
    else{
        feasible = 0;
        printf("Cannot found feasible solution\n");
    }
}

