/** @file
 * @brief Header file with internal function prototypes NCEPLIBS-g2c
 * library.
 *
 * @author Ed Hartnett @date 2021-11-08
 */

#ifndef _grib2_int_H
#define _grib2_int_H

#include <assert.h>
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <errno.h>
#if defined(WIN32)
#include <winsock2.h> /* ntohl() function for Windows. */
#else
#include <arpa/inet.h> /* ntohl() function for Unix/Mac. */
#endif

#ifdef PTHREADS
#include <pthread.h>
#endif

#include "grib2.h"

#define ALOG2 (0.69314718) /**< ln(2.0) */

/** Name of JPEG codec in Jasper. */
#define G2C_JASPER_JPEG_FORMAT_NAME "jpc"

/** Minimum acceptable value for max_bytes parameter of g2c_get_msg(). */
#define G2C_MIN_MAX_BYTES 16

/** GRIB magic header string. */
#define G2C_MAGIC_HEADER "GRIB"

/** Full length of magic header string (includes GRIB version byte). */
#define G2C_MAGIC_HEADER_LEN 8

/** Maximum number of messages in a file. */
#define G2C_MAX_MESSAGES 1024

#define BYTE 8  /**< Number of bits in a byte. */
#define WORD 32 /**< Number of bits in four bytes. */

#define ONE_BYTE 1    /**< One byte. */
#define TWO_BYTES 2   /**< Two bytes. */
#define FOUR_BYTES 4  /**< Four bytes. */
#define EIGHT_BYTES 8 /**< Eight bytes. */

/** Latitude */
#define LATITUDE "Latitude"

/** Longitude */
#define LONGITUDE "Longitude"

/* For GRIB1. */

/** Length of the product definition section (pds) in index file. */
#define G2C_INDEX1_PDS_VAL_LEN 27

/** Length of the grid definition section (gds) in index file. */
#define G2C_INDEX1_GDS_VAL_LEN 41

/** Length of the bitmap section (bms) in index file. */
#define G2C_INDEX1_BMS_VAL_LEN 5

/** Length of the binary data section (bds) in index file. */
#define G2C_INDEX1_BDS_VAL_LEN 10

/** Length of the pds 2 section in the index file. */
#define G2C_INDEX1_PDS_VAL2_LEN 27

/** Length of the pds 3 section in the index file. */
#define G2C_INDEX1_PDS_VAL3_LEN 27

/** Length of the gds in index file. */
#define G2C_INDEX1_GDS_VAL2_LEN 27

/* For thread-safety, use these macros. */
#ifdef PTHREADS

/** Define a mutex for thread-safety. */
#define MUTEX(m) pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;

/** Find an externally defined mutex for thread-safety. */
#define EXTERN_MUTEX(m) extern pthread_mutex_t m;

/** Define a macro to start a critical section for thread safety. */
#define MUTEX_LOCK(m)           \
    do                          \
    {                           \
        pthread_mutex_lock(&m); \
    } while (0)

/** Define a macro to end a critical section for thread safety. */
#define MUTEX_UNLOCK(m)           \
    do                            \
    {                             \
        pthread_mutex_unlock(&m); \
    } while (0)

#else /* not PTHREADS */

/** Pthreads not enabled, so do nothing. */
#define MUTEX(m)

/** Pthreads not enabled, so do nothing. */
#define EXTERN_MUTEX(m)

/** Pthreads not enabled, so do nothing. */
#define MUTEX_LOCK(m)

/** Pthreads not enabled, so do nothing. */
#define MUTEX_UNLOCK(m)

#endif /* PTHREADS */

/** Byte swap 64-bit ints. This converts big-endian 8-byte ints into
 * native endian 8-byte ints. */
#define ntoh64(y) (((uint64_t)ntohl(y)) << WORD | ntohl(y >> WORD))

/** Byte swap 64-bit ints. This converts native-endian 8-byte ints into
 * big-endian 8-byte ints. */
#define hton64(y) (((uint64_t)htonl(y)) << WORD | htonl(y >> WORD))

#define G2C_FILE_READ 0  /**< Read. */
#define G2C_FILE_WRITE 1 /**< Write. */

/** This is the information about each message. */
typedef struct g2c_message_info
{
    size_t msg_num;                       /**< Number of message in file (0-based). */
    size_t bytes_to_msg;                  /**< Number of bytes to skip in the file, to get to this message. */
    size_t bytes_in_msg;                  /**< Number of bytes in this message. */
    unsigned char discipline;             /**< Discipline from section 0. */
    int section1[G2C_SECTION1_ARRAY_LEN]; /**< Section 1 array. */
    int num_fields;                       /**< Number of fields in the message. */
    int num_local;                        /**< Number of local sections in the message. */
    size_t bytes_to_local;                /**< Number of bytes in the message before the (first) local section. */
    size_t bytes_to_bms;                  /**< Number of bytes in the message to the bitmap section. */
    size_t bytes_to_data;                 /**< Number of bytes in the message to the (first) data section. */
    int num_sections;                     /**< Number of sections in the file. */
    int *section_number;                  /**< Array (length num_sections) of section numbers. */
    size_t *section_offset;               /**< Array (length num_sections) of byte offsets from start of message to section. */
    int sec1_len;                         /**< Length of section 1. */
    short center;                         /**< Originating center. */
    short subcenter;                      /**< Originating subcenter. */
    unsigned char master_version;         /**< GRIB master tables version number. */
    unsigned char local_version;          /**< Version number of GRIB local tables used to augment Master Tables. */
    unsigned char sig_ref_time;           /**< Significance of reference time. */
    short year;                           /**< Year. */
    unsigned char month;                  /**< Month. */
    unsigned char day;                    /**< Day. */
    unsigned char hour;                   /**< Hour. */
    unsigned char minute;                 /**< Minute. */
    unsigned char second;                 /**< Second. */
    unsigned char status;                 /**< Production Status of Processed data in the GRIB message. */
    unsigned char type;                   /**< Type of processed data in this GRIB message. */
    struct g2c_section_info *sec;         /**< List of section metadata. */
    struct g2c_file_info *file;           /**< Pointer to containing file. */
    struct g2c_message_info *next;        /**< Pointer to next in list. */
} G2C_MESSAGE_INFO_T;

/** Information about a section 3 through 7 in a GRIB2 message. */
typedef struct g2c_section_info
{
    int sec_id;                    /**< ID of the section (0-based). */
    unsigned int sec_len;          /**< Length of the section (in bytes). */
    size_t bytes_to_sec;           /**< Number of bytes from start of message to this section. */
    unsigned char sec_num;         /**< Section number. */
    G2C_MESSAGE_INFO_T *msg;       /**< Pointer to contianing message. */
    void *sec_info;                /**< Pointer to struct specific for section 3, 4, 5, 6, or 7. */
    struct g2c_section_info *next; /**< Pointer to next in list. */
    struct g2c_section_info *prev; /**< Pointer to previous in list. */
    long long int *template;       /**< Grid, product, or data template. */
    int template_len;              /**< Number of entries in template. */
} G2C_SECTION_INFO_T;

/** Keep information about dimensions defined in section 3. */
typedef struct g2c_dim_info
{
    int dimid;                   /**< Dimension ID. */
    size_t len;                  /**< Length of dimension. */
    char name[G2C_MAX_NAME + 1]; /**< Name of dimension. */
    float *value;                /**< Array of dimension values. */
} G2C_DIM_INFO_T;

/** Information about [Section 3 GRID DEFINITION
 * SECTION](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect3.shtml). */
typedef struct g2c_section3_info
{
    unsigned char source_grid_def; /**< Source of grid definition (See Table 3.0). */
    unsigned int num_data_points;  /**< Number of data points. */
    unsigned char num_opt;         /**< Number of octets for optional list of numbers defining number of points. */
    unsigned char interp_list;     /**< Interpetation of list of numbers defining number of points (See Table 3.11). */
    unsigned short grid_def;       /**< Grid definition template number (= N) (See Table 3.1). */
    int *optional;                 /**< Optional list of numbers defining number of points. */
    G2C_DIM_INFO_T dim[2];         /**< Dimension information. */
} G2C_SECTION3_INFO_T;

/** Information about [Section 4 PRODUCT DEFINITION
 * SECTION](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect4.shtml). */
typedef struct g2c_section4_info
{
    int field_num;
    unsigned short num_coord; /**< Number of coordinate values after template. */
    unsigned short prod_def;  /**< Product definition template number (See Table 4.0). */
    int *optional;            /**< Optional list of numbers defining number of points. */
} G2C_SECTION4_INFO_T;

/** Information about [Section 5 DATA REPRESENTATION
 * SECTION](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect5.shtml). */
typedef struct g2c_section5_info
{
    /** Number of data points where one or more values are specified
     * in Section 7 when a bit map is present, total number of data
     * points when a bit map is absent. */
    unsigned int num_data_points;
    unsigned short data_def; /**< Data representation template number (See Table 5.0). */
} G2C_SECTION5_INFO_T;

/** Information about [Section 6 BIT-MAP
 * SECTION](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect6.shtml). */
typedef struct g2c_section6_info
{
    /** Bit map indicator. See [Table
     * 6.0](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table6-0.shtml). */
    unsigned char indicator;
} G2C_SECTION6_INFO_T;

/** This is the information about each open file. */
typedef struct g2c_file_info
{
    int g2cid;                   /**< ID of the file. */
    char path[G2C_MAX_NAME + 1]; /**< Path of the file. */
    FILE *f;                     /**< FILE pointer to open file. */
    size_t num_messages;         /**< Number of messages in the file. */
    G2C_MESSAGE_INFO_T *msg;     /**< Information about each message in the file. */
} G2C_FILE_INFO_T;

/** An entry in a GRIB2 code table. */
typedef struct g2c_entry
{
    struct g2c_entry *next;
    char desc[G2C_MAX_GRIB_DESC_LEN + 1];
    char status[G2C_MAX_GRIB_STATUS_LEN + 1];
    char code[G2C_MAX_GRIB_CODE_LEN + 1];
} G2C_CODE_ENTRY_T;

/** A GRIB2 code table. */
typedef struct g2c_code_table
{
    struct g2c_code_table *next;
    char title[G2C_MAX_GRIB_TITLE_LEN + 1];
    G2C_CODE_ENTRY_T *entry;
} G2C_CODE_TABLE_T;

/** An entry in the table of NOAA abbreviations for GRIB parameters. */
typedef struct g2c_param
{
    int g1num;                                /**< GRIB1 parameter number. */
    int g1ver;                                /**< GRIB1 table version. */
    int g2disc;                               /**< GRIB2 discipline. */
    int g2cat;                                /**< GRIB2 category number. */
    int g2num;                                /**< GRIB2 parameter number. */
    char abbrev[G2C_MAX_NOAA_ABBREV_LEN + 1]; /**< NOAA abbreviation for this parameter. */
} G2C_PARAM_T;

/**
 * Struct for GRIB template, returned by getgridtemplate().
 */
struct gtemplate
{
    /** The template type:
     * 3 Grid Defintion Template.
     * 4 Product Defintion Template.
     * 5 Data Representation Template. */
    g2int type;

    g2int num; /**< The template number. */

    /** Number of entries in the static part of the template. */
    g2int maplen;

    /** Number of octets of each entry in the static part of the
     * template. */
    g2int *map;

    /** Indicates whether or not the template needs to be extended. */
    g2int needext;

    /** Number of entries in the template extension. */
    g2int extlen;

    /** Number of octets of each entry in the extension part of the
     * template. */
    g2int *ext;
};

typedef struct gtemplate gtemplate; /**< Struct for GRIB template. */

/* Legacy support functions. */
double int_power(double x, g2int y);
void mkieee(float *a, g2int *rieee, g2int num);
void rdieee(g2int *rieee, float *a, g2int num);

/* Get the various templates. */
gtemplate *getdrstemplate(g2int number);
gtemplate *extdrstemplate(g2int number, g2int *list);
gtemplate *getpdstemplate(g2int number);
gtemplate *extpdstemplate(g2int number, g2int *list);
gtemplate *getgridtemplate(g2int number);
gtemplate *extgridtemplate(g2int number, g2int *list);

/* Level and date/time info. */
int g2c_get_datetime(int ipdtn, long long int *ipdtmpl, short year, unsigned char month, unsigned char day,
                     unsigned char hour, unsigned char minute, unsigned char second, char *tabbrev);
int g2c_get_level_desc(int ipdtn, long long int *ipdtmpl, char *level_desc);

/* Packing and unpacking data. */
void simpack(float *fld, g2int ndpts, g2int *idrstmpl,
             unsigned char *cpack, g2int *lcpack);
g2int simunpack(unsigned char *cpack, g2int *idrstmpl, g2int ndpts,
                float *fld);
void compack(float *fld, g2int ndpts, g2int idrsnum, g2int *idrstmpl,
             unsigned char *cpack, g2int *lcpack);
int comunpack(unsigned char *cpack, g2int lensec, g2int idrsnum,
              g2int *idrstmpl, g2int ndpts, float *fld);
void misspack(float *fld, g2int ndpts, g2int idrsnum, g2int *idrstmpl,
              unsigned char *cpack, g2int *lcpack);
void cmplxpack(float *fld, g2int ndpts, g2int idrsnum, g2int *idrstmpl,
               unsigned char *cpack, g2int *lcpack);
g2int getpoly(unsigned char *csec3, g2int *jj, g2int *kk, g2int *mm);
void specpack(float *fld, g2int ndpts, g2int JJ, g2int KK, g2int MM,
              g2int *idrstmpl, unsigned char *cpack, g2int *lcpack);
g2int specunpack(unsigned char *cpack, g2int *idrstmpl, g2int ndpts, g2int JJ,
                 g2int KK, g2int MM, float *fld);
g2int getdim(unsigned char *csec3, g2int *width, g2int *height, g2int *iscan);

int enc_png(unsigned char *data, g2int width, g2int height, g2int nbits,
            unsigned char *pngbuf);
int dec_png(unsigned char *pngbuf, g2int *width, g2int *height,
            unsigned char *cout);
void pngpack(float *fld, g2int width, g2int height, g2int *idrstmpl,
             unsigned char *cpack, g2int *lcpack);
g2int pngunpack(unsigned char *cpack, g2int len, g2int *idrstmpl, g2int ndpts,
                float *fld);
int enc_jpeg2000(unsigned char *cin, g2int width, g2int height, g2int nbits,
                 g2int ltype, g2int ratio, g2int retry, char *outjpc,
                 g2int jpclen);
int dec_jpeg2000(char *injpc, g2int bufsize, g2int *outfld);
void jpcpack(float *fld, g2int width, g2int height, g2int *idrstmpl,
             unsigned char *cpack, g2int *lcpack);
g2int jpcunpack(unsigned char *cpack, g2int len, g2int *idrstmpl, g2int ndpts,
                float *fld);

void aecpack(float *fld, g2int width, g2int height, g2int *idrstmpl,
             unsigned char *cpack, g2int *lcpack);
g2int aecunpack(unsigned char *cpack, g2int len, g2int *idrstmpl, g2int ndpts,
                float *fld);
int enc_aec(unsigned char *data, g2int ctemplen, g2int nbits, g2int flags,
            g2int block_size, g2int rsi, unsigned char *aecbuf, g2int *aecbuflen);
int dec_aec(unsigned char *cpack, g2int len, g2int nbits, g2int flags,
            g2int block_size, g2int rsi, unsigned char *cfld, g2int cfldlen);

/* Packing and unpacking bits. */
void gbit(unsigned char *in, g2int *iout, g2int iskip, g2int nbits);
void sbit(unsigned char *out, g2int *in, g2int iskip, g2int nbits);
void gbits(unsigned char *in, g2int *iout, g2int iskip, g2int nbits,
           g2int nskip, g2int n);
void sbits(unsigned char *out, g2int *in, g2int iskip, g2int nbits,
           g2int nskip, g2int n);

/* Deal with grib groups. */
int pack_gp(g2int *kfildo, g2int *ic, g2int *nxy,
            g2int *is523, g2int *minpk, g2int *inc, g2int *missp, g2int *misss,
            g2int *jmin, g2int *jmax, g2int *lbit, g2int *nov,
            g2int *ndg, g2int *lx, g2int *ibit, g2int *jbit, g2int *kbit,
            g2int *novref, g2int *lbitref, g2int *ier);

/* Check the message header and check for message termination. */
int g2c_check_msg(unsigned char *cgrib, g2int *lencurr, int verbose);

/* Basic file I/O. */
int g2c_file_io(FILE *f, int write, int g2ctype, void *var);
int g2c_file_io_byte(FILE *f, int write, char *var);
int g2c_file_io_ubyte(FILE *f, int write, unsigned char *var);
int g2c_file_io_short(FILE *f, int write, short *var);
int g2c_file_io_ushort(FILE *f, int write, unsigned short *var);
int g2c_file_io_int(FILE *f, int write, int *var);
int g2c_file_io_uint(FILE *f, int write, unsigned int *var);
int g2c_file_io_longlong(FILE *f, int write, long long *var);
int g2c_file_io_ulonglong(FILE *f, int write, unsigned long long *var);
int g2c_file_io_template(FILE *f, int rw_flag, int map, long long int *template_value);

/* Read and remember file, message, and section metadata. */
int g2c_add_file(const char *path, int mode, int *g2cid);
int add_msg(G2C_FILE_INFO_T *file, int msg_num, size_t bytes_to_msg, size_t bytes_in_msg,
            int read_file, G2C_MESSAGE_INFO_T **msg);
int add_section(FILE *f, G2C_MESSAGE_INFO_T *msg, int sec_id, unsigned int sec_len, size_t bytes_to_sec,
                unsigned char sec_num);
int g2c_rw_section1_metadata(FILE *f, int rw_flag, G2C_MESSAGE_INFO_T *msg);
int g2c_rw_section3_metadata(FILE *f, int rw_flag, G2C_SECTION_INFO_T *sec);
int g2c_rw_section4_metadata(FILE *f, int rw_flag, G2C_SECTION_INFO_T *sec);
int g2c_rw_section5_metadata(FILE *f, int rw_flag, G2C_SECTION_INFO_T *sec);
int g2c_rw_section6_metadata(FILE *f, int rw_flag, G2C_SECTION_INFO_T *sec);

/* Logging output. */
int g2c_log_section1(G2C_MESSAGE_INFO_T *msg);

/* Handle logging. */
#ifdef LOGGING

/* To log something... */
void g2_log(int severity, const char *fmt, ...);

/** Log a message to stdout. This is used for debugging the library. */
#define LOG(e) g2_log e

#else /* LOGGING */

/** Ignore logging to stdout. Library was not built with LOGGING=ON. */
#define LOG(e)

#endif /* LOGGING */

#endif /*  _grib2_int_H  */
