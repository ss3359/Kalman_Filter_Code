//
//  main.c
//  NATIONAL_PARKS_XCODE
//
//  Created by Owner on 1/14/26.
//
//
//#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>
//
//
//int Minimum(int arr[]){
//    int length=sizeof(arr)/sizeof(arr[0]);
//
//    int smallest=arr[0], index=0;
//    for(int i=0; i<length; i++){
//        if(arr[i]<=smallest){
//            index=i;
//        }
//    }
//    return index;
//}
//
//void ShortestLine(int n,int arrival_number){
//    int lines[n];
//   
//    
//    printf("Enter the number of people in %d lines \n", n);
//    
//    for (int i=0; i<n; i++){
//        scanf("%d \n",&lines[i]);
//    }
//
//    
//    while(arrival_number>0){
//        int shortest_line=Minimum(lines);
//        printf("Shortest at line %d with %d people \n", shortest_line,lines[shortest_line]);
//        lines[shortest_line]+=1;
//        arrival_number-=1;
//    }
//
//
//}
//
//int main(void){
//    int lines=3;
//    int arrival=4;
//
//    ShortestLine(lines, arrival);
//}
