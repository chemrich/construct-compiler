import argparse
import os
import threading
import time
import yaml
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from google import genai
from google.genai import types

# ---------------------------------------------------------------------------
# Rate limiter — prevents API throttling when running concurrent workers
# ---------------------------------------------------------------------------

class _RateLimiter:
    """Simple token-bucket rate limiter (thread-safe)."""

    def __init__(self, requests_per_minute: float):
        self._min_interval = 60.0 / requests_per_minute if requests_per_minute > 0 else 0
        self._lock = threading.Lock()
        self._last_request = 0.0

    def wait(self):
        """Block until it's safe to send the next request."""
        if self._min_interval <= 0:
            return
        with self._lock:
            now = time.monotonic()
            wait_time = self._last_request + self._min_interval - now
            if wait_time > 0:
                time.sleep(wait_time)
            self._last_request = time.monotonic()

PERSONAS = {
    "undergrad": {
        "system": "You are a skilled undergraduate student in a synthetic biology lab. You are curious and enthusiastic but sometimes mix up terminology or use informal names for parts (e.g., using 'GFP' instead of 'mEGFP', or confusing promoter types). Generate prompts for a genetic construct compiler that reflect this perspective.",
        "file": "prompt_corpus_undergrad.yaml"
    },
    "industry_vet": {
        "system": "You are a 20-year industry veteran who mostly makes AAV (Adeno-Associated Virus) vectors but is now dabbling with protein expression. You are used to AAV terminology and constraints. Generate prompts for a genetic construct compiler that reflect this perspective, often trying to translate AAV concepts to protein expression or asking for AAV-like constructs.",
        "file": "prompt_corpus_industry_vet.yaml"
    },
    "postdoc": {
        "system": "You are an overconfident postdoc with extensive experience in mammalian protein expression, but you have just started working in bacterial and yeast systems. You assume you know how everything works but might apply mammalian rules to bacteria/yeast inappropriately. Generate prompts for a genetic construct compiler that reflect this perspective.",
        "file": "prompt_corpus_postdoc.yaml"
    },
    "old_school": {
        "system": "You are an 'Old School' geneticist who has been cloning since the 1980s. You use outdated terminology, refer to parts by the name of the plasmid they were first discovered in, or use classic gene names instead of modern standard part names. Generate prompts for a genetic construct compiler that reflect this perspective.",
        "file": "prompt_corpus_old_school.yaml"
    },
    "minimalist": {
        "system": "You are an 'Extreme Minimalist' engineer obsessed with minimizing construct size. You actively try to omit parts that are usually required (like terminators or specific RBSs), assuming read-through or host machinery will compensate. Generate prompts for a genetic construct compiler that reflect this perspective.",
        "file": "prompt_corpus_minimalist.yaml"
    },
    "metabolic_engineer": {
        "system": "You are a 'Metabolic Engineer' building complex pathways. You want to express 4 to 6 enzymes in a specific order, often asking for 'operon style' expression in bacteria with balanced expression levels. Generate prompts for a genetic construct compiler that reflect this perspective.",
        "file": "prompt_corpus_metabolic_engineer.yaml"
    },
    "pure_logic": {
        "system": "You are a 'Pure Logic' synthetic biologist who thinks in terms of computer science and boolean logic. You describe constructs as gates, switches, and feedback loops rather than listing parts. Generate prompts for a genetic construct compiler that reflect this perspective.",
        "file": "prompt_corpus_pure_logic.yaml"
    },
    "repetition_tester": {
        "system": "You are a researcher who is 'Cloning Superstitious' and terrified of repetitive DNA sequences causing recombination. You ask for complex fusions or repeats but demand that no two sequences be identical, asking for sequence diversity or synonymous parts. Generate prompts for a genetic construct compiler that reflect this perspective.",
        "file": "prompt_corpus_repetition_tester.yaml"
    },
    "therapeutic": {
        "system": "You are a 'Therapeutic Developer' designing constructs for clinical use. You have strict 'must not have' constraints, such as forbidding viral sequences or antibiotic resistance genes. Generate prompts for a genetic construct compiler that reflect this perspective.",
        "file": "prompt_corpus_therapeutic.yaml"
    },
    "protein_engineer_bp": {
        "system": "You are a protein engineer working with *Bacillus* and *Pichia* expression systems. You are focused on optimizing protein yield, secretion, and stability. You often ask for specific secretion signals (like AmyE for Bacillus or alpha-factor for Pichia), codon optimization, or specific promoters suited for high-level expression in these hosts. Generate prompts for a genetic construct compiler that reflect this perspective.",
        "file": "prompt_corpus_protein_engineer_bp.yaml"
    }
}

def generate_prompts(persona, count, model_name="gemini-2.5-flash", existing_count=0, concurrency=1, rpm=0):
    api_key = os.environ.get("GEMINI_API_KEY") or os.environ.get("GOOGLE_API_KEY")
    if not api_key:
        raise ValueError("GOOGLE_API_KEY environment variable not set")
        
    client = genai.Client(api_key=api_key)
    
    persona_info = PERSONAS[persona]
    system_prompt = persona_info["system"]
    
    rate_limiter = _RateLimiter(rpm) if rpm > 0 else None
    
    chunk_size = 20
    chunks = (count + chunk_size - 1) // chunk_size
    
    all_prompts = []
    
    def _worker(i):
        start_idx = i * chunk_size
        current_chunk_size = min(chunk_size, count - start_idx)
        if current_chunk_size <= 0:
            return []
            
        print(f"  Generating chunk {i+1}/{chunks} ({current_chunk_size} prompts) for {persona}...")
        
        user_prompt = f"Generate {current_chunk_size} diverse prompts for testing a genetic construct compiler. Each prompt should describe a construct design request from your perspective. Return the results as a JSON list of objects, where each object has 'prompt' and 'expect' keys. The 'expect' key should be a dict containing expected properties like 'host' (e.g. 'e_coli', 'mammalian', 's_cerevisiae'), 'min_cistrons' (integer), and 'must_have_parts' (list of strings like 'PurificationTag', 'CleavageSite', 'SolubilityTag')."
        
        try:
            if rate_limiter:
                rate_limiter.wait()
                
            response = client.models.generate_content(
                model=model_name,
                contents=user_prompt,
                config=types.GenerateContentConfig(
                    response_mime_type="application/json",
                    system_instruction=system_prompt,
                ),
            )
            
            chunk_prompts = yaml.safe_load(response.text)
            if not isinstance(chunk_prompts, list):
                 print(f"Error: Expected list, got {type(chunk_prompts)}")
                 return []
            return chunk_prompts
        except Exception as e:
            print(f"Error calling or parsing Gemini: {e}")
            return []

    if concurrency <= 1:
        for i in range(chunks):
            all_prompts.extend(_worker(i))
    else:
        print(f"Running with {concurrency} parallel workers")
        with ThreadPoolExecutor(max_workers=concurrency) as executor:
            futures = {executor.submit(_worker, i): i for i in range(chunks)}
            for future in as_completed(futures):
                chunk_prompts = future.result()
                all_prompts.extend(chunk_prompts)
            
    processed_prompts = []
    for i, p in enumerate(all_prompts[:count]):
        p["id"] = f"{persona}_{existing_count + i}"
        p["category"] = f"persona_{persona}"
        p["difficulty"] = "medium"
        processed_prompts.append(p)
        
    return processed_prompts

def main():
    parser = argparse.ArgumentParser(description="Generate persona-based prompt corpus.")
    parser.add_argument("--persona", required=True, choices=PERSONAS.keys(), help="Persona to use")
    parser.add_argument("--count", type=int, default=500, help="Number of prompts to generate")
    parser.add_argument("--model", default="gemini-2.5-flash", help="Model to use")
    parser.add_argument("--concurrency", type=int, default=1, help="Number of parallel workers")
    parser.add_argument("--rpm", type=int, default=0, help="Rate limit in requests per minute")
    args = parser.parse_args()
    
    output_path = Path(__file__).parent / PERSONAS[args.persona]["file"]
    
    existing_prompts = []
    if output_path.exists():
        try:
            with open(output_path, 'r') as f:
                data = yaml.safe_load(f)
                if data and "prompts" in data:
                    existing_prompts = data["prompts"]
                    print(f"Found {len(existing_prompts)} existing prompts in {output_path}")
        except Exception as e:
            print(f"Error reading existing file: {e}")
            
    if len(existing_prompts) >= args.count:
        print(f"Already have {len(existing_prompts)} prompts. Nothing to do.")
        return
        
    needed = args.count - len(existing_prompts)
    print(f"Generating {needed} more prompts for persona: {args.persona} to reach {args.count}")
    
    new_prompts = generate_prompts(
        args.persona, 
        needed, 
        args.model, 
        existing_count=len(existing_prompts),
        concurrency=args.concurrency,
        rpm=args.rpm
    )
    
    all_prompts = existing_prompts + new_prompts
    output = {"prompts": all_prompts}
    
    with open(output_path, 'w') as f:
        yaml.dump(output, f, sort_keys=False)
        
    print(f"Successfully updated {output_path} with {len(new_prompts)} new prompts. Total: {len(all_prompts)}")

if __name__ == "__main__":
    main()
