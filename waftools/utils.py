import sys
from waflib.Configure import conf


@conf
def check_linux(conf):
    if not sys.platform.startswith('linux'):
        conf.fatal('Only Linux platform is supported at the moment')


def options(opt):
    opt.add_option(
        '--includepath',
        action='append',
        default=[],
        help='Additional include paths'
    )

    opt.add_option(
        '--cxxflags',
        action='append',
        default=[],
        help='Additional arguments for compiler'
    )

    opt.add_option(
        '--libpath',
        action='append',
        default=[],
        help='Additional paths to search for libraries'
    )


def configure(conf):
    conf.msg('Additional compiler flags', ' '.join(conf.options.cxxflags))
    conf.msg('Additional include paths', ','.join(conf.options.includepath))
    conf.msg('Additional paths to search libraries', ','.join(conf.options.libpath))

    conf.env.LIBPATH += conf.options.libpath
    conf.env.INCLUDES += conf.options.includepath
    conf.env.INCLUDES.append(conf.path.find_dir('src/libgcm').abspath())
    conf.env.INCLUDES.append(conf.path.find_dir('src').abspath())
    conf.env.CXXFLAGS += conf.options.cxxflags
    conf.env.LINKFLAGS += ['-lpthread', '-lrt']
    for x in conf.options.libpath:
        conf.env.LINKFLAGS += ['-Wl,-rpath,' + x]
