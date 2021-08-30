from pathlib import Path

from cobra.io import read_sbml_model
from cobra import Solution

from cobramod import add_reactions

# Load main model
backup = read_sbml_model(
    str(Path.cwd().joinpath("model", "PlantCoreMetabolism_v2_0_0.sbml"))
)
dir_data = Path.cwd().joinpath("data")
assert backup
assert dir_data.exists()


with open(file=Path.cwd().joinpath("model", "biomass.csv")) as f:
    f.readline()
    string = str()
    for line in f.readlines():
        identifier, coefficient = [string for string in line.split(",")][:2]
        coefficient = coefficient[1:]
        string += f"{coefficient} {identifier} + "
    string = f"New_biomass_tx,  New biomass for core model |{string[:-2]} -->"
model = backup.copy()

add_reactions(
    model=model, obj=string, directory=dir_data, show_imbalance=False
)
model.objective = "New_biomass_tx"
print(model.reactions.get_by_id("New_biomass_tx").reaction)

# Malate/Pyruvate transporter
add_reactions(
    model=model,
    obj="PYR_MAL_pc, Pyruvate/Malate transporter | "
    + "PYRUVATE_p + MAL_c <-> PYRUVATE_c + MAL_p",
    directory=dir_data,
)

assert model.reactions.get_by_id("PYR_MAL_pc")

# Initial import/export constraints - Defining autotrophic conditions
model.reactions.get_by_id("CO2_tx").bounds = (-1000, 1000)
model.reactions.get_by_id("H2O_tx").bounds = (-1000, 1000)
model.reactions.get_by_id("NH4_tx").bounds = (0.0, 0.0)
model.reactions.get_by_id("Pi_tx").bounds = (0, 1000)
model.reactions.get_by_id("SO4_tx").bounds = (0, 1000)
model.reactions.get_by_id("O2_tx").bounds = (-1000, 1000)
# Exported but not imported
model.reactions.get_by_id("Sucrose_tx").bounds = (-1000, 0)
model.reactions.get_by_id("GLC_tx").bounds = (-1000, 0)
# ATP bounds
model.reactions.get_by_id("ATPase_tx").bounds = (0, 1000)
# NTT is only active at night
model.reactions.get_by_id("ATP_ADP_Pi_pc").bounds = (0, 0)

c4_model = model.copy()
# Model below will be merge to c4_bundle
c4_bundle = model.copy()

# For Mesophyll cell
for dictlist in [
    c4_model.reactions,
    c4_model.metabolites,
    c4_model.genes,
    c4_model.groups,
]:
    for item in dictlist:
        item.id = f"M_{item.id}"

# For Bundle sheath cell
for dictlist in [
    c4_bundle.reactions,
    c4_bundle.metabolites,
    c4_bundle.genes,
    c4_bundle.groups,
]:
    for item in dictlist:
        item.id = f"B_{item.id}"

c4_model.merge(right=c4_bundle)

assert c4_model.reactions.get_by_id("M_CO2_tx")
assert c4_model.reactions.get_by_id("B_CO2_tx")
assert c4_model.reactions.get_by_id("M_ATPase_tx")
assert c4_model.reactions.get_by_id("B_ATPase_tx")

with open(file=Path.cwd().joinpath("model", "metabolites.txt")) as f:
    no_transport = [line[:-1].replace('"', "") for line in f.readlines()]

for metabolite in c4_model.metabolites:
    identifier = metabolite.id[2:]

    if identifier[: identifier.find("_")] in no_transport:
        continue
    elif f"MB_{identifier}" in [
        reaction.id for reaction in c4_model.reactions
    ]:
        continue

    elif "mc" not in metabolite.compartment and "c" in metabolite.compartment:
        add_reactions(
            model=c4_model,
            obj=f"MB_{identifier}, {identifier} bundle sheath/mesophyll cell "
            f"transporter | M_{identifier} <-> B_{identifier}",
            show_imbalance=False,
            directory=dir_data,
        )
        assert c4_model.reactions.get_by_id(f"MB_{identifier}")

# C4-specific constraints

# CONSTRAINT: No CO2 uptake in bundle sheath cells due to suberin layer
# in cell membranes
c4_model.reactions.get_by_id("B_CO2_tx").bounds = (0.0, 0.0)

# Force C4 cycle: Block Rubisco carboxylase/oxygenase in Mesophyll
c4_model.reactions.get_by_id(
    "M_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"
).bounds = (0, 0)

c4_model.reactions.get_by_id("M_RXN_961_p").bounds = (0, 0)

# Force NADP-ME decarboxylation pathway:
# Block all other decarboxylation reactions except NADP_ME in the plastid
c4_model.reactions.get_by_id("B_PEPCARBOXYKIN_RXN_c").bounds = (0, 0)
c4_model.reactions.get_by_id(
    "B_1_PERIOD_1_PERIOD_1_PERIOD_39_RXN_m"
).bounds = (0, 0)
c4_model.reactions.get_by_id("B_MALIC_NADP_RXN_c").bounds = (0, 0)

# Force NADP-ME decarboxylation pathways:
# make alternative decarboxylation routes irreversible
c4_model.reactions.get_by_id("B_CARBAMATE_KINASE_RXN_p").bounds = (
    0,
    1000,
)
c4_model.reactions.get_by_id("M_CARBAMATE_KINASE_RXN_p").bounds = (
    0,
    1000,
)

c4_model.reactions.get_by_id("B_ISOCITDEH_RXN_m").bounds = (0, 1000)
c4_model.reactions.get_by_id("M_ISOCITDEH_RXN_m").bounds = (0, 1000)

c4_model.reactions.get_by_id("B_ISOCITDEH_RXN_c").bounds = (0, 1000)
c4_model.reactions.get_by_id("M_ISOCITDEH_RXN_c").bounds = (0, 1000)

c4_model.reactions.get_by_id("B_ISOCITRATE_DEHYDROGENASE_NAD_RXN_m").bounds = (
    0,
    1000,
)

c4_model.reactions.get_by_id("M_ISOCITRATE_DEHYDROGENASE_NAD_RXN_m").bounds = (
    0,
    1000,
)


# Fix malate transport
c4_model.reactions.get_by_id("B_OAA_MAL_pc").bounds = (0, 1000)

c4_model.reactions.get_by_id("B_PYRUVATE_pc").bounds = (0, 0)

c4_model.reactions.get_by_id(
    "M_PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c"
).bounds = (0, 0)

# Constrains for light dependent maintenance costs
atp_b = c4_model.reactions.get_by_id("B_ATPase_tx")
photon_b = c4_model.reactions.get_by_id("B_Photon_tx")
atp_m = c4_model.reactions.get_by_id("M_ATPase_tx")
photon_m = c4_model.reactions.get_by_id("M_Photon_tx")

const_b = c4_model.problem.Constraint(
    (0.0049 * photon_b.flux_expression + 2.7852) - atp_b.flux_expression,
    lb=0,
    ub=0,
)
c4_model.add_cons_vars(const_b)

const_m = c4_model.problem.Constraint(
    (0.0049 * photon_m.flux_expression + 2.7852) - atp_m.flux_expression,
    lb=0,
    ub=0,
)
c4_model.add_cons_vars(const_m)

# ATP/NADPH 3:1 constraints
const = c4_model.problem.Constraint(
    c4_model.reactions.get_by_id("B_ATPase_tx").flux_expression
    - 3
    * (
        c4_model.reactions.get_by_id("B_NADPHoxc_tx").flux_expression
        + c4_model.reactions.get_by_id("B_NADPHoxp_tx").flux_expression
        + c4_model.reactions.get_by_id("B_NADPHoxm_tx").flux_expression
    ),
    lb=0,
    ub=0,
)
c4_model.add_cons_vars(const)

const = c4_model.problem.Constraint(
    c4_model.reactions.get_by_id("M_ATPase_tx").flux_expression
    - 3
    * (
        c4_model.reactions.get_by_id("M_NADPHoxc_tx").flux_expression
        + c4_model.reactions.get_by_id("M_NADPHoxp_tx").flux_expression
        + c4_model.reactions.get_by_id("M_NADPHoxm_tx").flux_expression
    ),
    lb=0,
    ub=0,
)
c4_model.add_cons_vars(const)


def c4_simulation(light, N, c4_model) -> Solution:
    with c4_model:
        # Light Uptale constrain
        B_Im_hnu = c4_model.reactions.get_by_id("B_Photon_tx")
        M_Im_hnu = c4_model.reactions.get_by_id("M_Photon_tx")

        # CONSTRAINT: Total Photon uptake limited to "light" µE
        const_hnu_sum = c4_model.problem.Constraint(
            B_Im_hnu.flux_expression + M_Im_hnu.flux_expression, lb=0, ub=light
        )
        c4_model.add_cons_vars(const_hnu_sum)

        # CONSTRAINT: Total Photon uptake by bundle sheath must be less or
        # equal than in mesophyll
        const_hnu_ratio = c4_model.problem.Constraint(
            M_Im_hnu.flux_expression - B_Im_hnu.flux_expression, lb=0, ub=light
        )
        c4_model.add_cons_vars(const_hnu_ratio)

        # CONSTRAINT : Total N uptake must not surpass defined upper bound
        bs_n = c4_model.reactions.get_by_id("B_Nitrate_tx")
        m_n = c4_model.reactions.get_by_id("M_Nitrate_tx")
        const_n_ratio = c4_model.problem.Constraint(
            bs_n.flux_expression + m_n.flux_expression, lb=0, ub=N
        )
        c4_model.add_cons_vars(const_n_ratio)

        return c4_bundle.optimize()
        # FBA


print
