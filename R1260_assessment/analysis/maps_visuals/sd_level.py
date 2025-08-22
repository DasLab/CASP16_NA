

def sd_level(session,maps,sdlevel):
    from numpy import float32
    from chimerax.core.commands import run
    for v in maps:
        m = v.full_matrix().astype(float32)
        print(m.std()*sdlevel)
        if m.max()>1e-9:
            run(session, f'volume #{v.id[0]} sdlevel {sdlevel}')

def register_command(session):
    from chimerax.core.commands import CmdDesc, register,FloatArg
    from chimerax.map import MapsArg
    desc = CmdDesc(required=[('maps', MapsArg),('sdlevel',FloatArg)],
                   synopsis='Set sdlevel is variance present')
    register('volume sd_level', desc, sd_level, logger=session.logger)

register_command(session)