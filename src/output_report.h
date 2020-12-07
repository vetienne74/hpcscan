#ifndef HPCSCAN_OUTPUT_REPORT_H_
#define HPCSCAN_OUTPUT_REPORT_H_

#include <string>

#include "type_def.h"

namespace hpcscan {

Rtn_code print_header_of_output_report(void) ;
Rtn_code print_end_of_output_report(Rtn_code) ;

void print_line1(void) ;
void print_line2(void) ;
void print_line3(void) ;
void print_line4(void) ;
void print_line5(void) ;
void print_blank(void) ;

void printInfo(Display_type, char*) ;
void printInfo(Display_type, const char*) ;
void printInfo(Display_type, char*, Myint) ;
void printInfo(Display_type, const char*, Myint) ;
void printInfo(Display_type, const char*, Myint, Myint) ;
void printInfo(Display_type, const char*, Myint64) ;
void printInfo(Display_type, char*, Myfloat) ;
void printInfo(Display_type, const char*, Myfloat32) ;
void printInfo(Display_type, const char*, Myfloat64) ;
void printInfo(Display_type, char*, char*) ;
void printInfo(Display_type, const char*, char*) ;
void printInfo(Display_type, const char*, string) ;
void printInfo(Display_type, const char*, const char*) ;

void printDebug(Debug_level, char*) ;
void printDebug(Debug_level, const char*) ;
void printDebug(Debug_level, char*, string) ;
void printDebug(Debug_level, char*, Myint) ;
void printDebug(Debug_level, char*, Myint, Myint) ;
void printDebug(Debug_level, char*, Myint64) ;
void printDebug(Debug_level, char*, Myfloat32) ;
void printDebug(Debug_level, char*, Myfloat64) ;
void printDebug(Debug_level, const char*, const char*) ;

void printWarning(string* text) ;
void printWarning(char* text) ;
void printWarning(char* text, Myint) ;
void printWarning(const char* text) ;

void printError(string* text) ;
void printError(char* text) ;
void printError(char* text, Myint) ;
void printError(char* text, string) ;
void printError(char* text, Myint, Myint) ;
void printError(char* text1, char* text2) ;
void printError(const char* text) ;

} // namespace hpcscan

#endif
