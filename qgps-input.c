#include "step.h"
#include "qgps-input.h"
#include <math.h>
#include <getopt.h>
#include <iniparser.h>
#include <stdio.h>

qgps_init_type_t qgps_init_type = QGPS_INIT_DELTA_K;
complex *qgps_init_data = NULL;
char *qgps_output_directory = NULL;
char *qgps_configuration_file = NULL;
int qgps_nx = 32;
int qgps_ny = 32;

const struct option options[] = {
        {"help", 0, 0, 'h'},
        {"configuration", 1, 0, 'f'},
        {"write-configuration", 1, 0, 'F'},
        {"initialize-type", 1, 0, 'i'},
        {"initialize-data", 1, 0, 'I'},
        {"timestep", 1, 0, 't'},
        {"domain-x", 1, 0, 'X'},
        {"domain-y", 1, 0, 'Y'},
        {"output-directory", 1, 0, 'D'},
        {0, 0, 0, 0}
};
char *option_string = "hf:F:I:i:t:X:Y:D:";
dictionary *config = NULL;
char *config_output = NULL;

int qgps_configure_broadcast();
int qgps_configure_master(int argc, char **argv);

int qgps_config_load(char *file);
int qgps_config_read();

char *qgps_option_get(const struct option *o);
int qgps_option_set(const struct option *o, char *value);
int qgps_option_read(const struct option *o);

const struct option * option_by_name(const char *name);
const struct option * option_by_flag(const char flag);
char *section();
char * sectioned_name(const struct option *o);
int number_of_options();

char next_option(char * const c, const int argc, char * const *argv, int * const index) {
        *c = getopt_long(argc, argv, option_string, options, index);
        *index = option_by_flag(*c) - options;

        return *c;
}

int qgps_configure(int argc, char **argv) {
        if (qgps_is_master_task)
                if (qgps_configure_master(argc, argv))
                        return 1;

        qgps_configure_broadcast();

        return 0;
}
int qgps_configure_master(int argc, char **argv) {
        int i = 0;
        for (char c = 0; next_option(&c, argc, argv, &i) != -1; ) {
                if (c == 'h') {
                        fprintf(stdout, "help output\n");
                        qgps_exit(EXIT_SUCCESS);
                        return 0;
                }
                else if (c == 'f')
                        qgps_config_load(optarg);
                else
                        qgps_option_set(&options[i], optarg);
        }
        return qgps_config_read();
}

int qgps_configure_broadcast() {
        MPI_Bcast(&qgps_init_type, sizeof(qgps_init_type_t),
                MPI_BYTE, qgps_master_task, QGPS_COMM_WORLD);

        MPI_Bcast(&qgps_time_step, 1,
                MPI_DOUBLE, qgps_master_task, QGPS_COMM_WORLD);
        MPI_Bcast(&qgps_nx, 1,
                MPI_DOUBLE, qgps_master_task, QGPS_COMM_WORLD);
        MPI_Bcast(&qgps_ny, 1,
                MPI_DOUBLE, qgps_master_task, QGPS_COMM_WORLD);

        int length = 0;
        if (qgps_is_master_task)
                length = strlen(qgps_output_directory);
        MPI_Bcast(&length, 1,
                MPI_INT, qgps_master_task, QGPS_COMM_WORLD);

        if (!qgps_is_master_task)
                qgps_output_directory = malloc((length + 1) * sizeof(char));
        MPI_Bcast(&qgps_output_directory, (length + 1) * sizeof(char),
                MPI_CHAR, qgps_master_task, QGPS_COMM_WORLD);

        return 0;
}

int qgps_config_load(char *file) {
        if (!config) {
                config = iniparser_load(file);
        }
        else {
                dictionary *new_conf = iniparser_load(file);
                if (!new_conf)
                        return 1;

                for (int i = 0; i < new_conf->n; i++)
                        qgps_option_set(option_by_name(new_conf->key[i]), new_conf->val[i]);

                iniparser_freedict(new_conf);
                return 0;
        }
}
const struct option *option_at(char *addr) {
      return (const struct option*)addr;
}
char qgps_configure_flag(char *name) {
        for (char *o = (char*)options; *o; o += sizeof(struct option))
                if (!strcmp(option_at(o)->name, name))
                        return option_at(o)->val;

        fprintf(stderr, "WARNING: Unrecognized option %s\n", name);
        return 0;
}
const struct option *option_by_name(const char *name) {
        const char *name2 = strstr(name, ":");
        if (name2)
                name = name2;

        for (char *o = (char*)options; *o; o += sizeof(struct option))
                if (!strcmp(option_at(o)->name, name))
                        return option_at(o);

        return NULL;
}
const struct option *option_by_flag(char flag) {
        for (char *o = (char*)options; *o; o += sizeof(struct option))
                if (option_at(o)->val == flag)
                        return option_at(o);

        return NULL;
}
int config_initialize() {
        config = dictionary_new(number_of_options());
        iniparser_setstring(config, section(), "");

        return 0;
}
int qgps_option_set(const struct option *o, char *value) {
        if (!o || !value)
                return 0;

        if (!config)
                config_initialize();

        iniparser_set(config, sectioned_name(o), value);
}
int qgps_config_read() {
        for (char *o = (char*)options; *o; o += sizeof(struct option))
                if (qgps_option_get(option_at(o)))
                        qgps_option_read(option_at(o));
}
char *qgps_option_get(const struct option *o) {
        if (!o)
                return 0;

        return iniparser_getstring(config, sectioned_name(o), NULL);
}
double qgps_option_getdouble(const struct option *o) {
        if (!o) {
                fprintf(stderr, "WARNING: getdouble() on NULL option; returning 0\n");
                return 0;
        }

        char *s = qgps_option_get(o), *end = s;
        double d = strtod(s, &end);
        if (*end || end == s) {
                fprintf(stderr, "WARNING: getdouble() failed on %s; returning NAN\n", o->name);
                return NAN;
        }

        return d;
}
int qgps_option_getint(const struct option *o) {
        if (!o) {
                fprintf(stderr, "WARNING: getint() on NULL option; returning 0\n");
                return 0;
        }

        char *s = qgps_option_get(o), *end = s;
        int i = strtol(s, &end, 10);
        if (*end || end == s) {
                fprintf(stderr, "WARNING: getint() failed on %s; returning 0\n", o->name);
                return 0;
        }

        return i;
}
int number_of_options() {
        static int n = -1;

        if (n < 0)
                for (n = 0; *((char*)&options[n]); n++);

        return n;
}
char *section() {
        return "configuration";
}
char *sectioned_name(const struct option *o) {
        static char **names = NULL;
        static int section_length = -1;

        if (!names)
                names = calloc(number_of_options(), sizeof(char*));
        if (section_length < 0)
                section_length = strlen(section());

        int i = o - options;
        if (!names[i]) {
                names[i] = malloc(section_length + strlen(o->name) + 2);

                strcpy(names[i], section());
                strcat(names[i], ":");
                strcat(names[i], o->name);
        }

        return names[i];
}

int qgps_option_read(const struct option *o) {
        if (!o)
                return 0;

        switch(o->val) {
        case 'h':
                fprintf(stdout, "help output\n");
                qgps_exit(EXIT_SUCCESS);
                break;
        case 'f':
                qgps_config_load(qgps_option_get(o));
                break;
        case 'F':
                config_output = qgps_option_get(o);
                break;
        case 'i':
                qgps_init_type = qgps_init_type_parse(qgps_option_get(o));
                break;
        case 'I':
                fprintf(stdout, "STUB: load initialization data from %s\n",
                        qgps_option_get(o));
                break;
        case 't':
                qgps_time_step = qgps_option_getdouble(o);
                break;
        case 'X':
                qgps_nx = qgps_option_getint(o);
                break;
        case 'Y':
                qgps_ny = qgps_option_getint(o);
                break;
        case 'd':
                qgps_output_directory = qgps_option_get(o);
                break;
        default:
                fprintf(stderr, "WARNING: Unimplemented option %s = %s", o->name, qgps_option_get(o));
                break;
        }

        return 0;
}

